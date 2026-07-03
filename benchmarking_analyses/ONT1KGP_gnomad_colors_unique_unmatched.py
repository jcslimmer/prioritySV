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

def extract_unmatched_sv_keys(fp_vcf):
    """
    Parses a Truvari false-positive VCF and returns the unmatched SV keys.
    """
    unmatched_keys = set()
    with pysam.VariantFile(fp_vcf) as vcf:
        for record in vcf:
            sv_key = get_sv_key(record)
            if sv_key:
                unmatched_keys.add(sv_key)
                
    print(f"Extracted {len(unmatched_keys)} unmatched SV keys from {fp_vcf}.")
    return unmatched_keys

def extract_intersected_svs(input_vcf, shared_unmatched_keys, output_vcf):
    """
    Parse VCF and extract SV records for unmatched keys.
    """
    records_written = 0
    with pysam.VariantFile(input_vcf) as vcf_in:
        with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:    
            for record in vcf_in:
                sv_key = get_sv_key(record)
                if sv_key and sv_key in shared_unmatched_keys:
                    vcf_out.write(record)
                    records_written += 1
    
    print(f"Wrote {records_written} unique unmatched SV records to {output_vcf}.")

if __name__ == "__main__":

    UNMATCHED_ONT1KGP_VCF = "./SV_database_matching/HG002_priority_merge_all_tp_comp_insertions_ONT1KGP/fp.vcf.gz"
    UNMATCHED_COLORS_VCF = "./SV_database_matching/HG002_priority_merge_all_tp_comp_insertions_gnomad/fp.vcf.gz"
    UNMATCHED_GNOMAD_VCF = "./SV_database_matching/HG002_priority_merge_all_tp_comp_insertions_colors/fp.vcf.gz"

    OUTPUT_VCF = "./SV_database_matching/HG002_priority_merge_all_unmatched_insertions.vcf"

    unmatched_keys_ONT1KGP = collect_unmatched_sv_keys(UNMATCHED_ONT1KGP_VCF)
    unmatched_keys_colors = collect_unmatched_sv_keys(UNMATCHED_COLORS_VCF)
    unmatched_keys_gnomad = collect_unmatched_sv_keys(UNMATCHED_GNOMAD_VCF)

    shared_unmatched_keys = unmatched_keys_ONT1KGP & unmatched_keys_colors & unmatched_keys_gnomad
    print(f"Found {len(shared_unmatched_keys)} structural variants that did not match any of the three databases.")

    extract_intersected_svs(UNMATCHED_ONT1KGP_VCF, shared_unmatched_keys, OUTPUT_VCF)