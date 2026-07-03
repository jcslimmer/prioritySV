import pysam

def get_sv_key(record):
    """
    Function to extract the comparison key from a VCF record.
    Returns (chrom, pos, abs(svlen)) or None if SVLEN is missing.
    """
    if 'SVLEN' in record.info:
        svlen = record.info['SVLEN']
        if isinstance(svlen, (list, tuple)):
            svlen = svlen[0]
        return (record.chrom, record.pos, abs(svlen))
    return None

def database_matchid_to_callset_sv_key_map(tp_comp_vcf):
    """
    Parses a tp-comp VCF to map the database SV MatchId (record in tp-base VCF)
    to the original callset's genomic coordinates.
    """
    mapping = {} # MatchId -> (chrom, pos, svlen)
    with pysam.VariantFile(tp_comp_vcf) as vcf:
        for record in vcf:
            if 'MatchId' in record.info:
                match_id = record.info['MatchId'][1]
                sv_key = get_sv_key(record)
                if sv_key:
                    mapping[match_id] = sv_key
    return mapping

def find_3way_shared_matches(manta_tp_base_vcf, gatksv_tp_base_vcf, priority_merge_tp_base_vcf,
    manta_tp_comp_vcf, gatksv_tp_comp_vcf, priority_merge_tp_comp_vcf, filters):
    """
    Finds the set of shared original calls (Manta, GATK-SV, priority_merge) that all matched 
    the same SV across each database. The matched database SVs must also all have an allele frequency >= 1%.

    Returns:
        dict: Key = (manta_sv_key, gatksv_sv_key, priority_merge_sv_key)
              Value = (manta_match_id, gatksv_match_id, priority_merge_match_id)
    """

    manta_matchid_sv_key_map = database_matchid_to_callset_sv_key_map(manta_tp_comp_vcf)
    gatksv_matchid_sv_key_map = database_matchid_to_callset_sv_key_map(gatksv_tp_comp_vcf)
    priority_merge_matchid_sv_key_map = database_matchid_to_callset_sv_key_map(priority_merge_tp_comp_vcf)

    # Database SVs matched to Manta callset 
    manta_base_matches = {} # database_sv_key -> manta_match_id
    with pysam.VariantFile(manta_tp_base_vcf) as vcf:

        has_af_header = 'AF' in vcf.header.info
        has_ac_header = 'AC' in vcf.header.info
        
        for record in vcf:
            # allele freq filter for gnomAD/CoLoRS, allele count for ONT1KGP
            if 'ac' in filters and has_ac_header:
                if record.info.get('AC', (0,))[0] < filters['ac']:
                    continue
            if 'af' in filters and has_af_header:
                if record.info.get('AF', (0,))[0] < filters['af']:
                    continue

            if 'MatchId' in record.info:
                sv_key = get_sv_key(record)
                if sv_key:
                    manta_base_matches[sv_key] = record.info['MatchId'][1]
    
    # Database SVs matched to GATK-SV callset 
    gatksv_base_matches = {}
    with pysam.VariantFile(gatksv_tp_base_vcf) as vcf:
        for record in vcf:
            if 'MatchId' in record.info:
                sv_key = get_sv_key(record)
                if sv_key and sv_key in manta_base_matches:
                    gatksv_base_matches[sv_key] = record.info['MatchId'][1]

    shared_base_sv_keys = []
    priority_merge_base_matches = {}
    with pysam.VariantFile(priority_merge_tp_base_vcf) as vcf:
        for record in vcf:
            if 'MatchId' in record.info:
                sv_key = get_sv_key(record)
                if sv_key and sv_key in gatksv_base_matches:
                    priority_merge_base_matches[sv_key] = record.info['MatchId'][1]
                    shared_base_sv_keys.append(sv_key)

    shared_sv_data = {} # (manta_comp_sv_key, gatksv_comp_sv_key, priority_merge_comp_sv_key) -> (manta_match_id, gatksv_match_id, priority_merge_match_id)

    for base_sv_key in shared_base_sv_keys:
        manta_match_id = manta_base_matches[base_sv_key]
        gatksv_match_id = gatksv_base_matches[base_sv_key]
        priority_merge_match_id = priority_merge_base_matches[base_sv_key]

        if (manta_match_id in manta_matchid_sv_key_map
            and gatksv_match_id in gatksv_matchid_sv_key_map
            and priority_merge_match_id in priority_merge_matchid_sv_key_map):
            
            shared_original_call_sv_key = (manta_matchid_sv_key_map[manta_match_id],
                                       gatksv_matchid_sv_key_map[gatksv_match_id],
                                       priority_merge_matchid_sv_key_map[priority_merge_match_id]
            )

            shared_sv_data[shared_original_call_sv_key] = (manta_match_id, gatksv_match_id, priority_merge_match_id)

    print(f"Found {len(shared_sv_data)} sets of Manta, GATK-SV and priority_merge calls matching the same >= 1% allele frequency database SV.")
    return shared_sv_data

def extract_shared_matches(tp_comp_vcf, final_comp_match_ids, output_vcf):
    """Extracts records from tp-comp VCF based on final set of MatchIds."""

    records_written = 0
    with pysam.VariantFile(tp_comp_vcf) as vcf_in:
        with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in.fetch():
                try:
                    if 'MatchId' in record.info:
                        match_id = record.info['MatchId'][1]
                        if match_id in final_comp_match_ids:
                            vcf_out.write(record)
                            records_written += 1
                except (KeyError, ValueError, TypeError, IndexError):
                    continue
    # should be equal for all output VCFs
    print(f"Wrote {records_written} records to {output_vcf}.")

if __name__ == "__main__":

    MANTA_TP_BASE_GNOMAD = "./HG002_manta_3way_shared_insertions_gnomad/tp-base.vcf.gz"
    GATKSV_TP_BASE_GNOMAD = "./HG002_gatksv_3way_shared_insertions_gnomad/tp-base.vcf.gz"
    PRIORITY_MERGE_TP_BASE_GNOMAD = "./HG002_priority_merge_3way_shared_insertions_gnomad/tp-base.vcf.gz"

    MANTA_TP_COMP_GNOMAD = "./HG002_manta_3way_shared_insertions_gnomad/tp-comp.vcf.gz"
    GATKSV_TP_COMP_GNOMAD = "./HG002_gatksv_3way_shared_insertions_gnomad/tp-comp.vcf.gz"
    PRIORITY_MERGE_TP_COMP_GNOMAD = "./HG002_priority_merge_3way_shared_insertions_gnomad/tp-comp.vcf.gz"

    MANTA_TP_BASE_ONT1KGP = "./HG002_manta_3way_shared_insertions_ONT1KGP/tp-base.vcf.gz"
    GATKSV_TP_BASE_ONT1KGP = "./HG002_gatksv_3way_shared_insertions_ONT1KGP/tp-base.vcf.gz"
    PRIORITY_MERGE_TP_BASE_ONT1KGP = "./HG002_priority_merge_3way_shared_insertions_ONT1KGP/tp-base.vcf.gz"

    MANTA_TP_COMP_ONT1KGP = "./HG002_manta_3way_shared_insertions_ONT1KGP/tp-comp.vcf.gz"
    GATKSV_TP_COMP_ONT1KGP = "./HG002_gatksv_3way_shared_insertions_ONT1KGP/tp-comp.vcf.gz"
    PRIORITY_MERGE_TP_COMP_ONT1KGP = "./HG002_priority_merge_3way_shared_insertions_ONT1KGP/tp-comp.vcf.gz"

    MANTA_TP_BASE_COLORS = "./HG002_manta_3way_shared_insertions_colors/tp-base.vcf.gz"
    GATKSV_TP_BASE_COLORS = "./HG002_gatksv_3way_shared_insertions_colors/tp-base.vcf.gz"
    PRIORITY_MERGE_TP_BASE_COLORS = "./HG002_priority_merge_3way_shared_insertions_colors/tp-base.vcf.gz"

    MANTA_TP_COMP_COLORS = "./HG002_manta_3way_shared_insertions_colors/tp-comp.vcf.gz"
    GATKSV_TP_COMP_COLORS = "./HG002_gatksv_3way_shared_insertions_colors/tp-comp.vcf.gz"
    PRIORITY_MERGE_TP_COMP_COLORS = "./HG002_priority_merge_3way_shared_insertions_colors/tp-comp.vcf.gz"

    shared_3way_matches_gnomad = find_3way_shared_matches(
        MANTA_TP_BASE_GNOMAD, GATKSV_TP_BASE_GNOMAD, PRIORITY_MERGE_TP_BASE_GNOMAD,
        MANTA_TP_COMP_GNOMAD, GATKSV_TP_COMP_GNOMAD, PRIORITY_MERGE_TP_COMP_GNOMAD,
        filters={'af: 0.01'}
        )
    
    shared_3way_matches_ont1kgp = find_3way_shared_matches(
        MANTA_TP_BASE_ONT1KGP, GATKSV_TP_BASE_ONT1KGP, PRIORITY_MERGE_TP_BASE_ONT1KGP,
        MANTA_TP_COMP_ONT1KGP, GATKSV_TP_COMP_ONT1KGP, PRIORITY_MERGE_TP_COMP_ONT1KGP,
        filters={'ac: 10'}
    )

    shared_3way_matches_colors = find_3way_shared_matches(
        MANTA_TP_BASE_COLORS, GATKSV_TP_BASE_COLORS, PRIORITY_MERGE_TP_BASE_COLORS,
        MANTA_TP_COMP_COLORS, GATKSV_TP_COMP_COLORS, PRIORITY_MERGE_TP_COMP_COLORS,
        filters={'af: 0.01'}
    )

    shared_database_matches = shared_3way_matches_gnomad.keys() & shared_3way_matches_ont1kgp.keys() & shared_3way_matches_colors.keys()
    print(f"Found {len(shared_database_matches)} sets of 3-way shared calls matching a >= 1% allele frequency SV in all 3 databases (gnomAD, ONT-1KGP and CoLoRS).")

    # Extract final 3-way shared and database intersected set of matched SVs from tp-comp VCFs
    manta_ids_gnomad = {shared_3way_matches_gnomad[k][0] for k in shared_database_matches}
    gatksv_ids_gnomad = {shared_3way_matches_gnomad[k][1] for k in shared_database_matches}
    priority_merge_ids_gnomad = {shared_3way_matches_gnomad[k][2] for k in shared_database_matches}
    
    extract_shared_matches(MANTA_TP_COMP_GNOMAD, manta_ids_gnomad, "HG002_manta_3way_shared_database_intersected_matches_insertions_gnomad.vcf")
    extract_shared_matches(GATKSV_TP_COMP_GNOMAD, gatksv_ids_gnomad, "HG002_gatksv_3way_shared_database_intersected_matches_insertions_gnomad.vcf")
    extract_shared_matches(PRIORITY_MERGE_TP_COMP_GNOMAD, priority_merge_ids_gnomad, "HG002_priority_merge_3way_shared_database_intersected_matches_insertions_gnomad.vcf")

    manta_ids_ont1kgp = {shared_3way_matches_ont1kgp[k][0] for k in shared_database_matches}
    gatksv_ids_ont1kgp = {shared_3way_matches_ont1kgp[k][1] for k in shared_database_matches}
    priority_merge_ids_ont1kgp = {shared_3way_matches_ont1kgp[k][2] for k in shared_database_matches}

    extract_shared_matches(MANTA_TP_COMP_ONT1KGP, manta_ids_ont1kgp, "HG002_manta_3way_shared_database_intersected_matches_insertions_ONT1KGP.vcf")
    extract_shared_matches(GATKSV_TP_COMP_ONT1KGP, gatksv_ids_ont1kgp, "HG002_gatksv_3way_shared_database_intersected_matches_insertions_ONT1KGP.vcf")
    extract_shared_matches(PRIORITY_MERGE_TP_COMP_ONT1KGP, priority_merge_ids_ont1kgp, "HG002_priority_merge_3way_shared_database_intersected_matches_insertions_ONT1KGP.vcf")

    manta_ids_colors = {shared_3way_matches_colors[k][0] for k in shared_database_matches}
    gatksv_ids_colors = {shared_3way_matches_colors[k][1] for k in shared_database_matches}
    priority_merge_ids_colors = {shared_3way_matches_colors[k][2] for k in shared_database_matches}

    extract_shared_matches(MANTA_TP_COMP_COLORS, manta_ids_colors, "HG002_manta_3way_shared_database_intersected_matches_insertions_colors.vcf")
    extract_shared_matches(GATKSV_TP_COMP_COLORS, gatksv_ids_colors, "HG002_gatksv_3way_shared_database_intersected_matches_insertions_colors.vcf")
    extract_shared_matches(PRIORITY_MERGE_TP_COMP_COLORS, priority_merge_ids_colors, "HG002_priority_merge_3way_shared_database_intersected_matches_insertions_colors.vcf")
