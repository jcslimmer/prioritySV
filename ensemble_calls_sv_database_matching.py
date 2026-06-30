import argparse
import pysam
import re

def infer_input_vcf(header):
    """
    Infer the SV database VCF by searching the headers for a database-specific identifier.
    """

    header_strings = str(header)
    sv_database=None

    if 'gnomAD_SV' in header_strings:
        sv_database = "gnomAD"
    elif 'pbsv' in header_strings:
        sv_database = "CoLoRS"
    elif 'Sniffles' in header_strings:
        sv_database = "ONT_1KGP"
    
    return sv_database

def create_new_header(input_vcf, new_info_fields=None):
    """ 
    Create a new, modified VCF header.

    Args:
        input_vcf (str): Path to the original VCF file.
        new_info_fields (list of dicts, optional): A list of new INFO fields to add.
            Each dict should be {'id': 'ID', 'number': 'Num', 'type': 'Type', 'desc': 'Desc'}.
    """
    with pysam.VariantFile(input_vcf, 'r') as vcf_reader:
        original_header = vcf_reader.header
        new_header = original_header.copy()

        # Add new database INFO fields
        if new_info_fields:
            for field in new_info_fields:
                new_header.info.add(field['id'], field['number'], field['type'], field['desc'])
    
    return new_header

def copy_pysam_record(record, vcf_out):
    """
    Creates a new pysam record and copies all data from an input record.
    """
    new_record = vcf_out.new_record()
    new_record.chrom = record.chrom
    new_record.pos = record.pos
    new_record.id = record.id
    new_record.qual = record.qual

    for f in record.filter:
        new_record.filter.add(f)
    
    if record.alleles is None:
        svtype = record.info.get('SVTYPE')
        symbolic_alt = f"<{svtype}>"
        ref_base = record.ref if record.ref is not None else 'N'
        new_record.alleles = (ref_base, symbolic_alt)
    else:
        new_record.alleles = record.alleles
    
    for key, value in record.info.items():
        new_record.info[key] = value
    
    for sample_name in record.samples:
        for key, value in record.samples[sample_name].items():
            new_record.samples[sample_name][key] = value

    return new_record

def annotate_combine_vcfs(tp_base_vcf, tp_comp_vcf, fp_vcf, output_vcf, sv_database):
    """
    Finds SVs shared between two VCFs by chrom, pos, svlen, and collects their MatchIds.

    Args:
        tp_base_vcf (str): Path to callset VCF.
        tp_comp_vcf (str): Path to comparison VCF output by Truvari.

    """   
    # Build lookup table from callset VCF
    allele_freq_dict = {}

    with pysam.VariantFile(tp_base_vcf) as vcf_in:
        for base_record in vcf_in.fetch():
            if 'MatchId' in base_record.info and 'AC' in base_record.info:
                base_vcf_id = base_record.info['MatchId']
                comp_vcf_id = base_vcf_id[1]

                allele_count = base_record.info['AC'][0]
                if 'ONT_1KGP' in sv_database:
                    allele_freq = allele_count / 900
                else:
                    allele_freq = base_record.info['AF'][0]
                locus = f"{base_record.chrom}:{base_record.pos}-{base_record.stop}"
                database_ref_seq = base_record.ref
                database_alt_seq = base_record.alts[0]
                allele_freq_dict[comp_vcf_id] = allele_freq, locus, database_ref_seq, database_alt_seq

    print(f"{len(allele_freq_dict)} SVs from callset matched to {sv_database} database SVs.")

    new_info_fields = [
        {'id': f'{sv_database}_AF', 'number': 'A', 'type': 'Float', 'desc': f'Allele frequency in {sv_database}'},
        {'id': f'{sv_database}_COORDS', 'number': '1', 'type': 'String', 'desc': f'Genomic coordinates of {sv_database} matched SV (1-based)'},
        {'id': f'{sv_database}_REF', 'number': '1', 'type': 'String', 'desc': f'Reference sequence for {sv_database} matched SV'},
        {'id': f'{sv_database}_ALT', 'number': '1', 'type': 'String', 'desc': f'Alternate sequence for {sv_database} matched SV'},
        {'id': f'{sv_database}_PctSizeSimilarity', 'number': '1', 'type': 'Float', 'desc': f'Percent size similarity between this variant and its closest match in {sv_database}'},
        {'id': f'{sv_database}_PctRecOverlap', 'number': '1', 'type': 'Float', 'desc': f'Percent reciprocal overlap with its closest match in {sv_database}'},
        {'id': f'{sv_database}_SizeDiff', 'number': '1', 'type': 'Integer', 'desc': f'Size difference from the {sv_database} matched SV'},
        {'id': f'{sv_database}_StartDistance', 'number': '1', 'type': 'Integer', 'desc': f'StartDistance difference from the {sv_database} matched SV'},
        {'id': f'{sv_database}_EndDistance', 'number': '1', 'type': 'Float', 'desc': f'EndDistance difference from the {sv_database} matched SV'}
    ]
    new_header = create_new_header(tp_comp_vcf, new_info_fields=new_info_fields)

    # Extract original calls from TP-comp VCF and recreate record with new INFO fields
    with pysam.VariantFile(output_vcf, 'w', header=new_header) as vcf_out:
        with pysam.VariantFile(tp_comp_vcf) as vcf_in:
            for comp_record in vcf_in.fetch():
                new_comp_record = copy_pysam_record(comp_record, vcf_out)
                new_record_alleles = None
                svlen = comp_record.info.get('SVLEN')[0]
                svtype = comp_record.info.get('SVTYPE')

                if 'MatchId' in new_comp_record.info:
                    comp_match_id_tuple = new_comp_record.info['MatchId']
                    comp_match_id = comp_match_id_tuple[1]
                    if comp_match_id in allele_freq_dict:
                        
                        allele_freq, locus, database_ref_seq, database_alt_seq = allele_freq_dict[comp_match_id]

                        new_comp_record.info[f'{sv_database}_AF'] = allele_freq
                        new_comp_record.info[f'{sv_database}_COORDS'] = locus
                        new_comp_record.info[f'{sv_database}_REF'] = database_ref_seq
                        new_comp_record.info[f'{sv_database}_ALT'] = database_alt_seq

                new_comp_record.info[f'{sv_database}_PctSizeSimilarity'] = new_comp_record.info.get('PctSizeSimilarity')
                new_comp_record.info[f'{sv_database}_PctRecOverlap'] = new_comp_record.info.get('PctRecOverlap')
                new_comp_record.info[f'{sv_database}_SizeDiff'] = new_comp_record.info.get('SizeDiff')
                new_comp_record.info[f'{sv_database}_StartDistance'] = new_comp_record.info.get('StartDistance')
                new_comp_record.info[f'{sv_database}_EndDistance'] = new_comp_record.info.get('EndDistance')

                if 'INS' in svtype:
                    new_comp_record.stop = new_comp_record.pos
                if 'DUP' in svtype:
                    new_comp_record.stop = new_comp_record.pos + svlen
                if 'DEL' in svtype:
                    new_comp_record.stop = new_comp_record.pos + abs(svlen)
                if 'INV' in svtype:
                    new_comp_record.stop = new_comp_record.pos + svlen   
                
                vcf_out.write(new_comp_record)
    # Extract records with no database matches and recreate record to handle new header INFO fields
        with pysam.VariantFile(fp_vcf) as vcf_in:
            for fp_record in vcf_in.fetch():
                new_fp_record = copy_pysam_record(fp_record, vcf_out)
                svlen = new_fp_record.info.get('SVLEN')[0]
                svtype = new_fp_record.info.get('SVTYPE')

                new_fp_record.info[f'{sv_database}_AF'] = None
                new_fp_record.info[f'{sv_database}_COORDS'] = None
                new_fp_record.info[f'{sv_database}_REF'] = None
                new_fp_record.info[f'{sv_database}_ALT'] = None
                new_fp_record.info[f'{sv_database}_PctSizeSimilarity'] = None
                new_fp_record.info[f'{sv_database}_PctRecOverlap'] = None
                new_fp_record.info[f'{sv_database}_SizeDiff'] = None
                new_fp_record.info[f'{sv_database}_StartDistance'] = None
                new_fp_record.info[f'{sv_database}_EndDistance'] = None
                
                if 'INS' in svtype:
                    new_fp_record.stop = new_fp_record.pos
                if 'DUP' in svtype:
                    new_fp_record.stop = new_fp_record.pos + svlen
                if 'DEL' in svtype:
                    new_fp_record.stop = new_fp_record.pos + abs(svlen)
                if 'INV' in svtype:
                    new_fp_record.stop = new_fp_record.pos + svlen
 
                vcf_out.write(new_fp_record)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotates a SV VCF callset with allele frequencies using database matched SV calls from truvari bench."
    )

    parser.add_argument(
        "-b", "--base_vcf",
        metavar="PATH",
        required=True,
        help="Path to a base VCF containing allele count (AC) info (tp-base.vcf.gz) from bcftools +fill-tags."
    )
    parser.add_argument(
        "-c", "--comp_vcf",
        metavar="PATH",
        required=True,
        help="Path to the comp VCF with original SVs passing matching thresholds (tp-comp.vcf.gz)."
    )
    parser.add_argument(
        "-f", "--fp_vcf",
        metavar="PATH",
        required=True,
        help="Path to the fp VCF with SVs not passing matching thresholds (fp.vcf.gz)."
    )
    parser.add_argument(
        "-o", "--output_vcf",
        metavar="PATH",
        required=True,
        help="Path to allele frequency annotated output VCF file."
    )

    args = parser.parse_args()
    
    with pysam.VariantFile(args.base_vcf, 'r') as vcf_reader:
        sv_database = infer_input_vcf(vcf_reader.header)
    
    annotate_combine_vcfs(args.base_vcf, args.comp_vcf, args.fp_vcf, args.output_vcf, sv_database)