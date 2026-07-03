#!/usr/bin/env python3

import argparse
import sys
import pysam
import csv
import os
from typing import Optional, Dict

# Standard hg38 contig lengths
HG38_CONTIGS = [
    ('chr1', 248956422), ('chr2', 242193529), ('chr3', 198295559),
    ('chr4', 190214555), ('chr5', 181538259), ('chr6', 170805979),
    ('chr7', 159345973), ('chr8', 145138636), ('chr9', 138394717),
    ('chr10', 133797422), ('chr11', 135086622), ('chr12', 133275309),
    ('chr13', 114364328), ('chr14', 107043718), ('chr15', 101991189),
    ('chr16', 90338345), ('chr17', 83257441), ('chr18', 80373285),
    ('chr19', 58617616), ('chr20', 64444167), ('chr21', 46709983),
    ('chr22', 50818468), ('chrX', 156040895), ('chrY', 57227415),
]

CHROM_ORDER = {f'chr{i}': i for i in range(1, 23)}
CHROM_ORDER['chrX'] = 23
CHROM_ORDER['chrY'] = 24

def get_sort_key(chrom: str) -> int:
    """Returns a sort key for karyotypic chromosome sorting."""
    return CHROM_ORDER.get(chrom, 99)

def infer_input_format(input_file: str) -> Optional[str]:
    """
    Infer the file format by checking the first few lines of the input file.
    """
    try:
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')

                if 'matchclips' in line or line.startswith("#CHR"):
                    return 'matchclip'
                if 'MEI_Family' in line:
                    return 'insertion_txt'
                if 'DEL.START' in line:
                    return 'deletion_txt'
                   
            return None
    except Exception:
        return None

def create_new_header(sample_name: str) -> pysam.VariantHeader:
    """ 
    Create minimal VCF header for Truvari benchmarking.
    """
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")

    for contig, length in HG38_CONTIGS:
        header.add_line(f"##contig=<ID={contig},length={length}>")

    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">')
    header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    header.add_line('##FORMAT=<ID=GT,Number=G,Type=String,Description="Genotype">')
    
    header.add_sample(sample_name)
    return header

def parse_info_string(info_str: str) -> Dict[str, str]:
    """Parses a VCF INFO string into a dictionary."""
    info_dict = {}
    if not info_str or info_str == '.': 
        return info_dict
    fields = info_str.split(';')
    for field in fields:
        if '=' in field:
            key, value = field.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[field] = "True"
    return info_dict

def process_insertion_deletion_txt(input_file: str, output_vcf: str, new_header: pysam.VariantHeader, file_format: str, sample_name: str) -> bool:
    """
    Process PredictedDeletion.txt and MEIs.txt files output by Scramble.
    """
    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    parse_error_count = 0
    svlen_filter_count = 0
    total_pass_records_written = 0

    records_buffer = []
    auto_sex_chromosomes = {contig for contig, length in HG38_CONTIGS}

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        header_row = next(reader)
        header_map = {name: i for i, name in enumerate(header_row)}

        for i, row in enumerate(reader):
            vcf_records_processed += 1
            try:
                variant_data = None
                chrom = ""
                pos = 0
                end = 0
                svlen = 0
                svtype = ""
                ref = 'N'
                
                if file_format == 'insertion_txt':
                    chrom_pos = row[header_map['Insertion']]
                    chrom, pos_str = chrom_pos.split(':')
                    
                    if chrom not in auto_sex_chromosomes:
                        non_auto_sex_chromosome_variants += 1
                        continue
                        
                    pos = int(pos_str)
                    end = pos
                    ins_start_len = int(row[header_map['Start_In_MEI']])
                    ins_stop_len = int(row[header_map['Stop_In_MEI']])
                    svlen = ins_stop_len - ins_start_len
                    svtype = "INS"
                    
                elif file_format == 'deletion_txt':
                    chrom = row[header_map['CONTIG']]
                    
                    if chrom not in auto_sex_chromosomes:
                        non_auto_sex_chromosome_variants += 1
                        continue

                    pos = int(row[header_map['DEL.START']])
                    end = int(row[header_map['DEL.END']])
                    svlen = int(row[header_map['DEL.LENGTH']])
                    ref = row[header_map['REF.ANCHOR.BASE']]
                    svtype = "DEL"
                    svlen = -abs(svlen)

                if abs(svlen) < 50:
                    svlen_filter_count += 1
                    continue

                variant_data = {
                    'chrom': chrom, 'pos': pos, 'end': end, 
                    'svtype': svtype, 'svlen': svlen, 'ref': ref
                }
                records_buffer.append(variant_data)

            except (ValueError, IndexError, KeyError):
                parse_error_count += 1
                continue

    # Sort records
    records_buffer.sort(key=lambda r: (get_sort_key(r['chrom']), r['pos']))

    # Write records
    vcf_writer = pysam.VariantFile(output_vcf, 'w', header=new_header)

    for variant in records_buffer:
        sv_record = vcf_writer.new_record()
        sv_record.contig = variant['chrom']
        sv_record.pos = variant['pos']
        sv_record.stop = variant['end']
        sv_record.ref = variant['ref']
        sv_record.alts = (f"<{variant['svtype']}>",)
        
        sv_record.info['SVTYPE'] = variant['svtype']
        sv_record.info['SVLEN'] = variant['svlen']
        
        sv_record.filter.add('PASS')
        sv_record.samples[sample_name]['GT'] = (0, 1)
        sv_record.samples[sample_name].phased = False
        
        vcf_writer.write(sv_record)
        total_pass_records_written += 1

    vcf_writer.close()

    print(f"{vcf_records_processed} records parsed from input file.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} records skipped because SVLEN < 50bp.")
    print(f"{parse_error_count} records skipped due to parsing/formatting errors.")
    print(f"{total_pass_records_written} sorted VCF records written to {output_vcf}.")

    return True

def process_matchclip_data(input_file: str, output_vcf: str, new_header: pysam.VariantHeader, file_format: str, sample_name: str) -> bool:
    """
    Processes matchclip data output by MOPline in the {sample}.matchclip.out file.
    """
    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    parse_error_count = 0
    svlen_filter_count = 0
    total_pass_records_written = 0

    records_buffer = []
    auto_sex_chromosomes = {contig for contig, length in HG38_CONTIGS}

    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            row = line.split('\t')
            vcf_records_processed += 1
            
            try:
                variant_data = None
                chrom = ""

                chrom, pos_str, end_str, svtype, svlen_str = row[:5]
                
                if chrom not in auto_sex_chromosomes:
                    non_auto_sex_chromosome_variants += 1
                    continue
                
                pos, end = int(pos_str), int(end_str)
                svlen = -int(svlen_str) if svtype in ("DUP", "DEL") else int(svlen_str)
                
                variant_data = { 'chrom': chrom, 'pos': pos, 'end': end, 'svtype': svtype, 'svlen': svlen }

                if variant_data:
                    records_buffer.append(variant_data)

            except (ValueError, IndexError):
                parse_error_count += 1
                continue

    # Sort records
    records_buffer.sort(key=lambda r: (get_sort_key(r['chrom']), r['pos']))

    # Write records
    vcf_writer = pysam.VariantFile(output_vcf, 'w', header=new_header)

    for variant in records_buffer:
        sv_record = vcf_writer.new_record()
        sv_record.contig = variant['chrom']
        sv_record.pos = variant['pos']
        sv_record.stop = variant['end']
        sv_record.ref = 'N'
        sv_record.alts = (f"<{variant['svtype']}>",)
        
        sv_record.info['SVTYPE'] = variant['svtype']
        sv_record.info['SVLEN'] = variant['svlen']
        
        sv_record.filter.add('PASS')
        sv_record.samples[sample_name]['GT'] = (0, 1)
        sv_record.samples[sample_name].phased = False
        
        vcf_writer.write(sv_record)
        total_pass_records_written += 1

    vcf_writer.close()

    print(f"{vcf_records_processed} records parsed from input file.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not 1-22, X or Y. Records excluded.")
    print(f"{parse_error_count} records skipped due to parsing/formatting errors.")
    print(f"{total_pass_records_written} sorted VCF records written to {output_vcf}.")

    return True

def main():
    parser = argparse.ArgumentParser(
        description="Converts MOPline raw MATCHCLIP output and (PredictedDeletions.txt, MEIs.txt - Scramble) to a minimal VCF compatible with Truvari."
    )
    parser.add_argument(
        "input_file", 
        help="Path to the input file."
    )
    parser.add_argument(
        "output_vcf", 
        help="Path for the output VCF file."
    )
    parser.add_argument(
        "--sample-name", 
        default="SAMPLE", 
        help="Sample name to use for the output VCF."
    )
    
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: The input file was not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    file_format = infer_input_format(args.input_file)
    if not file_format:
        print(f"Error: Could not determine the format of '{args.input_file}'.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Inferred file format: '{file_format}'", flush=True)

    new_header = create_new_header(args.sample_name)

    try:
        if file_format in ["insertion_txt", "deletion_txt"]:
            process_insertion_deletion_txt(args.input_file, args.output_vcf, new_header, file_format, args.sample_name)
        elif file_format  == "matchclip":
            process_matchclip_data(args.input_file, args.output_vcf, new_header, file_format, args.sample_name)
        else:
             print(f"Error: Inferred format '{file_format}' is not supported.", file=sys.stderr)
             sys.exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()