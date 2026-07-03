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

def collect_match_ids(manta_vcf_path, gatksv_vcf_path, priority_merge_vcf_path):
    """
    Finds SVs shared between two VCFs by chrom, pos, svlen, and collects their MatchIds.

    Args:
        manta_vcf_path (str): Path to the Manta tp-base VCF.
        gatksv_vcf_path (str): Path to the GATK-SV tp-base VCF.
        priority_merge_vcf_path (str): Path to the priority_merge tp-base VCF.

    Returns:
        tuple: A tuple containing three sets (manta_ids, gatksv_ids, priority_merge_ids).
    """
    print("Finding detected SVs shared and extracting caller/pipeline respective tp-comp MatchIds")
    
    # Build lookup table from the Manta VCF
    manta_sv = {}
    with pysam.VariantFile(manta_vcf_path) as vcf_in:
        for record in vcf_in.fetch():
            if 'MatchId' in record.info:
                sv_key = get_sv_key(record)
                if sv_key:
                    manta_sv[sv_key] = record.info['MatchId']

    print(f"Found {len(manta_sv)} unique variants in Manta VCF.")

    # Build lookup table from GATK-SV VCF
    gatksv_sv = {}
    with pysam.VariantFile(gatksv_vcf_path) as vcf_in:
        for record in vcf_in.fetch():
            if 'MatchId' in record.info:
                sv_key = get_sv_key(record)
                if sv_key:
                    gatksv_sv[sv_key] = record.info['MatchId']
    
    print(f"Found {len(gatksv_sv)} unique variants in GATK-SV VCF.")

    # Store respective unique MatchIds
    manta_ids = set()
    gatksv_ids = set()
    priority_merge_ids = set()
    seen_manta_ids = set()

    # 3-way intersection to find shared Manta, GATK-SV and priority_merge detected truth set SVs
    with pysam.VariantFile(priority_merge_vcf_path) as vcf_in:
        for record in vcf_in.fetch():
            if 'MatchId' in record.info:
                sv_key = get_sv_key(record)
                if sv_key:
                    # check if SV was also found by Manta and GATK-SV
                    if sv_key in manta_sv and sv_key in gatksv_sv:
                        manta_match_tuple = manta_sv[sv_key]
                        gatksv_match_tuple = gatksv_sv[sv_key]
                        priority_merge_match_tuple = record.info['MatchId']
                        # check if Manta MatchId already seen
                        manta_comp_id = manta_match_tuple[1]
                        if manta_comp_id not in seen_manta_ids:
                        # add the callset-specific MatchId to the sets
                            manta_ids.add(manta_match_tuple[1])
                            gatksv_ids.add(gatksv_match_tuple[1])
                            priority_merge_ids.add(priority_merge_match_tuple[1])
                            seen_manta_ids.add(manta_comp_id)
                        else:
                            continue

    # Check unique MatchIds are equal between callsets
    print(f"Discovered {len(manta_ids)} unique Manta MatchIds from shared variants.")
    print(f"Discovered {len(gatksv_ids)} unique GATK-SV MatchIds from shared variants.")
    print(f"Discovered {len(priority_merge_ids)} unique priority_merge MatchIds in 3-way intersection.")

    return manta_ids, gatksv_ids, priority_merge_ids

def extract_records_by_match_id(input_vcf_path, match_ids, output_vcf_path):
    """
    Reads an input VCF and writes records to a new VCF if their MatchId is in a given set.

    Args:
        input_vcf_path (str): Path to the source VCF to filter.
        match_ids (set): A set of MatchIds to look for.
        output_vcf_path (str): Path to write the filtered output VCF file.
    """
    
    records_written = 0

    with pysam.VariantFile(input_vcf_path) as vcf_in:
        with pysam.VariantFile(output_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in.fetch():
                try:
                    record_match_id = record.info['MatchId'][1]
                    # write record if MatchId is in the set of ids
                    if record_match_id in match_ids:
                        vcf_out.write(record)
                        records_written += 1
                except (KeyError, ValueError, IndexError):
                    # skip records without MatchIds or malformed MatchIds
                    continue
                    
    print(f"Wrote {records_written} records to: {output_vcf_path}")

if __name__ == "__main__":

    # TP base files output by Truvari bench
    MANTA_TP_BASE_VCF =  "./HG002_manta_dupins_refdist500_size70/HG002_manta_dupins_refdist500_size70_tp-base.vcf"
    GATKSV_TP_BASE_VCF = "./HG002_gatksv_dupins_refdist500_size70/HG002_gatksv_dupins_refdist500_size70_tp-base.vcf"
    PRIORITY_MERGE_TP_BASE_VCF = "./HG002_priority_merge_dupins_refdist500_size70/tp-base.vcf.gz"

    # TP comp VCF files output by Truvari bench
    MANTA_TP_COMP_VCF = "./HG002_manta_dupins_refdist500_size70/tp-comp_repeat_annotated.vcf.gz"
    GATKSV_TP_COMP_VCF = "./HG002_gatksv_dupins_refdist500_size70/tp-comp_repeat_annotated.vcf.gz"
    PRIORITY_MERGE_TP_COMP_VCF = "./HG002_priority_merge_dupins_refdist500_size70/tp-comp_repeat_annotated.vcf.gz"

    # Output files for truth set calls shared by Manta, GATK-SV and priority_merge
    MANTA_OUTPUT_VCF = "HG002_manta_3way_shared_insertions.vcf"
    GATKSV_OUTPUT_VCF = "HG002_gatksv_3way_shared_insertions.vcf"
    PRIORITY_MERGE_OUTPUT_VCF = "HG002_priority_merge_3way_shared_insertions.vcf"
    
    manta_ids, gatksv_ids, priority_merge_ids = collect_match_ids(MANTA_TP_BASE_VCF, GATKSV_TP_BASE_VCF, PRIORITY_MERGE_TP_BASE_VCF)

    extract_records_by_match_id(
        input_vcf_path=MANTA_TP_COMP_VCF,
        match_ids=manta_ids,
        output_vcf_path=MANTA_OUTPUT_VCF
    )

    extract_records_by_match_id(
        input_vcf_path=GATKSV_TP_COMP_VCF,
        match_ids=gatksv_ids,
        output_vcf_path=GATKSV_OUTPUT_VCF
    )

    extract_records_by_match_id(
        input_vcf_path=PRIORITY_MERGE_TP_COMP_VCF,
        match_ids=priority_merge_ids,
        output_vcf_path=PRIORITY_MERGE_OUTPUT_VCF
    )