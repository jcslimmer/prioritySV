#!/usr/bin/env python3

import argparse
import sys
import pysam
import numpy as np
import re
import os
from typing import Dict

HG38_CONTIGS = [
    ('chr1', 248956422), ('chr2', 242193529), ('chr3', 198295559),
    ('chr4', 190214555), ('chr5', 181538259), ('chr6', 170805979),
    ('chr7', 159345973), ('chr8', 145138636), ('chr9', 138394717),
    ('chr10', 133797422), ('chr11', 135086622), ('chr12', 133275309),
    ('chr13', 114364328), ('chr14', 107043718), ('chr15', 101991189),
    ('chr16', 90338345), ('chr17', 83257441), ('chr18', 80373285),
    ('chr19', 58617616), ('chr20', 64444167), ('chr21', 46709983),
    ('chr22', 50818468), ('chrX', 156040895), ('chrY', 57227415)
]

def get_sort_key(chrom):
    """Function for sorting chromosomes (1-22, X=23, Y=24)."""
    sort_chrom = chrom.replace("chr", "")
    if sort_chrom.isdigit():
        return int(sort_chrom)
    if sort_chrom == "X":
        return 23
    if sort_chrom == "Y":
        return 24
    
    return 99

def parse_info_string(info_str):
    """Parses VCF INFO string from a line stripped and tab-split VCF record into a dictionary."""
    info_dict = {}
    if not info_str or info_str == '.':
        return info_dict
    
    fields = info_str.split(';')
    for field in fields:
        if '=' in field:
            key, value = field.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[field] = True
    return info_dict

def create_minimal_vcf_header(sample_name):
    """
    Create minimal header for MOPline SV caller headerless VCF.
    """
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")

    for contig, length in HG38_CONTIGS:
        header.add_line(f"##contig=<ID={contig},length={length}>")

    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of the structural variant">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">')

    header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    header.add_line('##FORMAT=<ID=GT,Number=G,Type=String,Description="Genotype">')
    header.add_sample(sample_name)
    return header

def process_mopline_headerless_vcf(input_vcf, output_sv_vcf, sample_name):
    """
    Processes MOPline SV caller headerless VCF. 
    """
    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    svlen_filter_count = 0
    no_svlen_svtype_count = 0
    total_records_written = 0
    
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    
    minimal_header = create_minimal_vcf_header(sample_name)

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, 'w', header=minimal_header)
    vcf_reader = open(input_vcf, 'r')

    for i, record in enumerate(vcf_reader):
        if record.startswith('#'):
            continue
        record = record.strip()
        if not record:
            continue
        
        vcf_records_processed += 1

        vcf_parts = record.strip().split('\t')

        chrom = vcf_parts[0]
        pos_str = vcf_parts[1]
        filter_str = vcf_parts[6]
        info_str = vcf_parts[7]
        

        if chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue

        pos = int(pos_str)
        info_dict = parse_info_string(info_str)
        
        svtype = info_dict.get('SVTYPE')
        svlen = int(info_dict.get('SVLEN'))
        
        if not svtype or svlen is None:
            no_svlen_svtype_count += 1
            continue

        if abs(svlen) < 50:
            svlen_filter_count += 1
            continue

        if svtype in ('DEL', 'DUP', 'INV'):
            end = pos + abs(svlen)
        else: 
            end = pos

        sv_record = vcf_writer_sv.new_record()
        sv_record.contig = chrom
        sv_record.pos = pos
        sv_record.stop = end
        sv_record.ref = 'N'
        sv_record.alts = (f'<{svtype}>',)
                
        sv_record.info['SVTYPE'] = svtype
        sv_record.info['SVLEN'] = svlen
                                    
        sv_record.samples[sample_name]['GT'] = (0, 1)
        sv_record.samples[sample_name].phased = False
        
        vcf_writer_sv.write(sv_record)
        total_records_written += 1
    
    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} records skipped because SVLEN < 50bp.")
    print(f"{no_svlen_svtype_count} records skipped due to no SVLEN or SVTYPE.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

def infer_input_vcf(header, first_record):
    """
    Infer the SV caller VCF by searching the headers for a caller-specific identifier
    """

    header_strings = str(header)

    if 'INSurVeyorCommand' in header_strings:
        return 'insurveyor'
    elif 'DYSGU' in header_strings:
        return 'dysgu'
    elif 'SurVIndel2Command' in header_strings:
        return 'survindel2'
    elif 'GenerateSVCandidates' in header_strings:
        return 'manta'
    elif 'gcSV' in header_strings:
        return 'gcsv'
    elif 'LUMPY probability' in header_strings:
        return 'lumpy_smoove'
    elif 'MELT' in header_strings:
        return 'melt'
    elif 'WHAM' in header_strings:
        return 'wham'
    elif 'CNVnator' in header_strings:
        return 'cnvnator'
    elif 'SURVIVOR' in header_strings:
        return 'parliament'
    elif 'MOPline-1.7' in header_strings:
        return 'mopline'
    elif 'CPX' in header_strings:
        return 'gatksv'
    elif 'SVMETHOD' in header.info and first_record is not None:
        svmethod = first_record.info.get('SVMETHOD')
        if svmethod is not None and 'DELLY' in svmethod:
            return 'delly'
    else:
        return None

def create_new_header(input_vcf, sv_caller, new_sample_name=None):
    """ 
    Create new VCF header.
    - Change SVLEN Number to . for compatibility with other SV caller VCFs.
    - Change sample name for compatibility with sample name present in other SV caller VCFs
    - Changes specific SurVeyor headers for compatibility with pysam.
    """

    vcf_reader = pysam.VariantFile(input_vcf, 'r')
    original_header = vcf_reader.header
    new_header = pysam.VariantHeader()

    ids_to_replace = {'SVLEN', 'GQ', 'SVINSSEQ', 'GT', 'END'}
  
    for line in str(original_header).splitlines():
        match = re.search(r'##[A-Z]+=<ID=([^,>]+)', line)
        if match and match.group(1) in ids_to_replace:
            continue
        if line.startswith('#CHROM'):
            continue
        new_header.add_line(line)

    new_header.info.add('SVLEN', '.', 'Integer', 'Difference in length between REF and ALT alleles')
    new_header.info.add('END', '1', 'Integer', 'End position of the variant described in this record.')
    new_header.formats.add('GQ', '1', 'Integer', 'Genotype Quality')
    new_header.formats.add('GT', 'G', 'String', 'Genotype')
    new_header.info.add('SVINSSEQ', '.', 'String', 'Sequence of insertion')
    new_header.info.add('UNMERGED_SVTYPE', '1', 'String', 'Type of structural variant before Truvari bench --dup-to-ins')
    new_header.info.add('UNRESOLVED_DUPINS', '0', 'Flag', 'Sequence of DUP/INS is not fully assembled. ALT sequence represents concatenation of sequence assembled at both breakpoints.')

    if new_sample_name:
        new_header.add_sample(new_sample_name)
    else:
    # Check for duplicated sample names and deduplicate if necessary
        for sample_name in original_header.samples:
            sample_name_dup = sample_name.split('_')
            if len(sample_name_dup) == 2 and sample_name_dup[0] == sample_name_dup[1]:
                standardised_sample_name = sample_name_dup[0]
            else:
                standardised_sample_name = sample_name

            new_header.add_sample(standardised_sample_name)

    vcf_reader.close()
    return new_header

def process_mopline_vcf(vcf_reader, output_sv_vcf, new_header, sv_caller_name):
    
    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    svlen_filter_count = 0
    symbolic_ins_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)
    
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value
        svtype = record.info.get('SVTYPE')

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                # Round genotype qualities to nearest integer to match VCF header format
                if key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value


        write_sv_record = True

        if 'INS' in svtype or sv_record.alts[0] == "<INS>":
            sv_record.info['SVLEN'] = svlen
            symbolic_ins_count += 1

        elif 'DUP' in svtype or sv_record.alts[0] == "<DUP>":
            sv_record.info['SVLEN'] = svlen
            sv_record.stop = sv_record.pos + abs(svlen)
            symbolic_dup_count += 1

        elif 'DEL' in svtype or sv_record.alts[0] == "<DEL>":
            sv_record.info['SVLEN'] = -svlen
            sv_record.stop = sv_record.pos + abs(svlen)
            symbolic_del_count += 1
        else:
            sv_record.info['SVLEN'] = svlen
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{inv_svtype_count} VCF records SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_breakdancer_vcf(input_vcf, output_sv_vcf, sample_name):
    """
    Process malformed Breakdancer VCF from Parliament2 output for compatibility with Truvari bench/collapse.
    """
    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    no_svlen_svtype_count = 0
    svlen_filter_count = 0
    total_records_written = 0
    
    records_buffer = []
    minimal_header = create_minimal_vcf_header(sample_name)
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    vcf_reader = open(input_vcf, 'r')
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, 'w', header=minimal_header)

    for i, record in enumerate(vcf_reader):
        if record.startswith('#'):
            continue
        record = record.strip()
        if not record:
            continue

        vcf_records_processed += 1

        vcf_parts = record.split('\t')
        chrom = vcf_parts[0]
        pos_str = vcf_parts[1]
        info_str = vcf_parts[7]

        if chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue

        pos = int(pos_str)
        info_dict = parse_info_string(info_str)
        svtype = info_dict.get('SVTYPE')
        
        if svtype not in ('DEL', 'INV'):
            no_svlen_svtype_count += 1
            continue
        
        end = int(info_dict.get('END', 0))
        
        if 'INV' in svtype:
            svlen = (end - pos)
        if 'DEL' in svtype:
            svlen = -(end - pos)
        
        if abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        
        # Buffer data for sorting before writing
        variant_data = {
            'chrom': chrom,
            'pos': pos, 
            'end': end,
            'svtype': svtype,
            'svlen': svlen,
        }
        records_buffer.append(variant_data)

    # Sort records
    records_buffer.sort(key=lambda r: (get_sort_key(r['chrom']), r['pos']))

    # Create and write records
    for variant in records_buffer:
        sv_record = vcf_writer_sv.new_record()
        sv_record.contig = variant['chrom']
        
        sv_record.pos = variant['pos'] 
        sv_record.stop = variant['end']
        sv_record.ref = 'N'
        sv_record.alts = (f"<{variant['svtype']}>",)
               
        sv_record.samples[sample_name]['GT'] = (0, 1)
        sv_record.samples[sample_name].phased = False

        sv_record.info['SVTYPE'] = variant['svtype']
        sv_record.info['SVLEN'] = variant['svlen']
    
        vcf_writer_sv.write(sv_record)
        total_records_written += 1
    
    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} records skipped because SVLEN < 50bp.")
    print(f"{no_svlen_svtype_count} records skipped due to no SVLEN or SVTYPE.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

def process_parliament_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    symbolic_ins_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    new_header.formats.add('SP', '1', 'String', 'Supporting SV callers for this sample')

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    info_fields = ['SVTYPE', 'END', 'SUPP', 'CALLERS']
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = int(np.round(value[0])) if isinstance(value := record.info.get('AVGLEN'), tuple) else int(np.round(value))

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        # Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if 'TRA' in svtype or 'BND' in svtype:
            bndtra_record = vcf_writer_bndtra.new_record()
            bndtra_record.chrom = record.chrom
            bndtra_record.pos = record.pos
            bndtra_record.id = f"{sv_caller_name}_{record.id}"
            bndtra_record.alleles = record.alleles
            bndtra_record.qual = record.qual
            bndtra_record.filter.clear()
            for f in record.filter:
                bndtra_record.filter.add(f)
      
            for key in info_fields:
                if key in record.info:
                    bndtra_record.info[key] = record.info[key]

            for i, sample in enumerate(record.samples):
                sample_data = record.samples[i]
                for key, value in sample_data.items():
                    if key == 'GT':
                        sv_record.samples[i]['GT'] = value

            vcf_writer_bndtra.write(bndtra_record)
            passed_bndtra_records_written += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.qual = record.qual
        sv_record.filter.clear()
        for f in record.filter:
            sv_record.filter.add(f)

        for key in info_fields:
            if key in record.info:
                sv_record.info[key] = record.info[key]

        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                if key == 'GT':
                    sv_record.samples[i]['GT'] = value

        write_sv_record = True

        if 'INS' in svtype or sv_record.alts[0] == "<INS>":
            sv_record.info['SVLEN'] = svlen
            symbolic_ins_count += 1

        elif 'DUP' in svtype or sv_record.alts[0] == "<DUP>":
            sv_record.info['SVLEN'] = svlen
            sv_record.stop = sv_record.pos + abs(svlen)
            symbolic_dup_count += 1

        elif 'DEL' in svtype or sv_record.alts[0] == "<DEL>":
            sv_record.info['SVLEN'] = -svlen
            sv_record.stop = sv_record.pos + abs(svlen)
            symbolic_del_count += 1
        else:
            sv_record.info['SVLEN'] = svlen
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()
    
    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF record SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{inv_svtype_count} VCF records SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

def process_wham_cnvnator_vcf(vcf_reader, output_sv_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    svlen_filter_count = 0
    total_records_written = 0
    
    for contig, length in HG38_CONTIGS:
        new_header.add_line(f"##contig=<ID={contig},length={length}>")

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                if value is not None:
                    sv_record.samples[i][key] = value
        
        svtype = sv_record.info['SVTYPE']
        sv_record.alts = (f"<{svtype}>",)

        if svtype in ('DEL','INV', 'DUP') and svlen is not None:
            sv_record.stop = sv_record.pos + abs(svlen)
        if 'INS' in svtype:
            sv_record.stop = sv_record.pos

        vcf_writer_sv.write(sv_record)
        total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_melt_vcf(vcf_reader, output_sv_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    svlen_filter_count = 0
    total_records_written = 0

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)
    
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                # Round genotype qualities to nearest integer to match VCF header format
                if key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value

        sv_record.alts = ("<INS>",)
        sv_record.info['SVTYPE'] = "INS"
        vcf_writer_sv.write(sv_record)
        total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_lumpy_smoove_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    non_symbolic_dup_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    non_symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        # Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if 'TRA' in svtype or 'BND' in svtype:
            bndtra_record = vcf_writer_bndtra.new_record()
            bndtra_record.chrom = record.chrom
            bndtra_record.pos = record.pos
            bndtra_record.id = f"{sv_caller_name}_{record.id}"
            bndtra_record.alleles = record.alleles
            bndtra_record.info.update(record.info)

            for i, sample in enumerate(record.samples):
                sample_data = record.samples[i]
                for key, value in sample_data.items():
                    if key == 'GQ' and value is not None:
                        bndtra_record.samples[i][key] = int(np.round(float(value)))
                    elif value is not None:
                        bndtra_record.samples[i][key] = value
      
            vcf_writer_bndtra.write(bndtra_record)
            passed_bndtra_records_written += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                # Round genotype qualities to nearest integer to match VCF header format
                if key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value

        write_sv_record = True

        if 'INS' in svtype:
            sv_record.stop = sv_record.pos
            if sv_record.alts[0] == "<INS>":
                symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1

        elif 'DUP' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            sv_record.info['UNMERGED_SVTYPE'] = "DUP"
            if sv_record.alts[0] == '<DUP>':
                sv_record.info['UNRESOLVED_DUPINS'] = True
                symbolic_dup_count += 1
            else:
                non_symbolic_dup_count += 1

        elif 'DEL' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == "<DEL>":
                symbolic_del_count += 1
            else:
                non_symbolic_del_count += 1
        else:
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_delly_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    non_symbolic_dup_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    non_symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value
        end = record.info.get('END')
        pos = record.pos

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        if svtype is None:
            continue
        # Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if 'TRA' in svtype or 'BND' in svtype:

            bnd_pattern = re.compile(r'[\[\]](.+?):')
            bnd_chr_match = bnd_pattern.search(record.alts[0])
            bnd_chr = bnd_chr_match.group(1)
            if bnd_chr not in auto_sex_chromosomes:
                continue
            
            bndtra_record = vcf_writer_bndtra.new_record()
            bndtra_record.chrom = record.chrom
            bndtra_record.pos = record.pos
            bndtra_record.id = f"{sv_caller_name}_{record.id}"
            bndtra_record.alleles = record.alleles
            bndtra_record.info.update(record.info)

            for i, sample in enumerate(record.samples):
                sample_data = record.samples[i]
                for key, value in sample_data.items():
                    if key == 'GQ' and value is not None:
                        bndtra_record.samples[i][key] = int(np.round(float(value)))
                    elif value is not None:
                        bndtra_record.samples[i][key] = value
      
            vcf_writer_bndtra.write(bndtra_record)
            passed_bndtra_records_written += 1
            continue
        # Calculate SVLEN from END coordinate if it exists otherwise skip records with no SVLEN
        if svlen is None:
            if 'INS' not in svtype:
                if end is not None and isinstance(end, (list, tuple)):
                    end = end[0]
                if end is None:
                    end = record.stop
                if end is not None:
                    svlen = -(end-pos) if 'DEL' in svtype else end-pos
                    record.info['SVLEN'] = svlen
            else:
                print(f"Skipping {svtype} at {record.chrom}:{record.pos} - No END found")
                continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                # Round genotype qualities to nearest integer to match VCF header format
                if key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value

        write_sv_record = True

        if 'INS' in svtype:
            sv_record.stop = sv_record.pos
            if sv_record.alts[0] == "<INS>":
                symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1

        elif 'DUP' in svtype:
            sv_record.info['UNMERGED_SVTYPE'] = "DUP"
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == '<DUP>':
                sv_record.info['UNRESOLVED_DUPINS'] = True
                symbolic_dup_count += 1
            else:
                non_symbolic_dup_count += 1

        elif 'DEL' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == "<DEL>":
                symbolic_del_count += 1
            else:
                non_symbolic_del_count += 1
        else:
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True


def process_gcsv_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    non_symbolic_dup_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    non_symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value
        end = record.info.get('END')
        pos = record.pos

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        if svtype is None:
            continue
        # Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if 'TRA' in svtype or 'BND' in svtype:

            bnd_pattern = re.compile(r'[\[\]](.+?):')
            bnd_chr_match = bnd_pattern.search(record.alts[0])
            bnd_chr = bnd_chr_match.group(1)
            if bnd_chr not in auto_sex_chromosomes:
                continue
            
            bndtra_record = vcf_writer_bndtra.new_record()
            bndtra_record.chrom = record.chrom
            bndtra_record.pos = record.pos
            bndtra_record.id = f"{sv_caller_name}_{record.id}"
            bndtra_record.alleles = record.alleles
            bndtra_record.info.update(record.info)

            for i, sample in enumerate(record.samples):
                sample_data = record.samples[i]
                for key, value in sample_data.items():
                    if key == 'GQ' and value is not None:
                        bndtra_record.samples[i][key] = int(np.round(float(value)))
                    elif value is not None:
                        bndtra_record.samples[i][key] = value
      
            vcf_writer_bndtra.write(bndtra_record)
            passed_bndtra_records_written += 1
            continue

        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                # Round genotype qualities to nearest integer to match VCF header format
                if key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value

        write_sv_record = True

        if 'INS' in svtype:
            sv_record.stop = sv_record.pos
            if sv_record.alts[0] == "<INS>":
                symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1

        elif 'DUP' in svtype:
            sv_record.info['UNMERGED_SVTYPE'] = "DUP"
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == '<DUP>':
                sv_record.info['UNRESOLVED_DUPINS'] = True
                symbolic_dup_count += 1
            else:
                non_symbolic_dup_count += 1

        elif 'DEL' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == "<DEL>":
                symbolic_del_count += 1
            else:
                non_symbolic_del_count += 1
        else:
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_gatksv_vcf(vcf_reader, output_sv_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    cnv_svtype_count = 0
    svlen_filter_count = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    non_symbolic_dup_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    non_symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    info_fields = ['SVTYPE', 'END', 'SVLEN']

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value
        end = record.info.get('END')
        pos = record.pos

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        if svtype is None:
            continue
        if 'CNV' in svtype:
            cnv_svtype_count += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.qual = record.qual
        for f in record.filter:
            sv_record.filter.add(f)

        for key in info_fields:
            if key in record.info:
                sv_record.info[key] = record.info[key]

        # Iterate over samples
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                if key == 'GT':
                    sv_record.samples[i]['GT'] = value

        write_sv_record = True

        if 'INS' in svtype:
            sv_record.stop = sv_record.pos
            if sv_record.alts[0] == "<INS>":
                symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1

        elif 'DUP' in svtype:
            sv_record.info['UNMERGED_SVTYPE'] = "DUP"
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == '<DUP>':
                sv_record.info['UNRESOLVED_DUPINS'] = True
                symbolic_dup_count += 1
            else:
                non_symbolic_dup_count += 1

        elif 'DEL' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            if sv_record.alts[0] == "<DEL>":
                symbolic_del_count += 1
            else:
                non_symbolic_del_count += 1
        else:
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{cnv_svtype_count} VCF records SVTYPE is CNV. Records excluded.")
    print(f"{svlen_filter_count} VCF records with SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_insurveyor_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    symbolic_ins_unassembled_surveyor_seq_count = 0
    symbolic_ins_assembled_seq_count = 0
    symbolic_ins_unassembled_seq_count = 0
    total_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        ins_seq = record.info.get('SVINSSEQ')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        ### add counter
        if svtype is None:
            continue
        # Skip INS with no SVLEN and Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if svlen is None:
            if 'INS' in svtype:
                continue
            if 'TRA' in svtype or 'BND' in svtype:
                if 'dysgu' in sv_caller_name:
                    continue
                alt_allele = str(record.alts[0])

                bnd_pattern = re.compile(r'[\[\]](.+?):')
                bnd_chr_match = bnd_pattern.search(alt_allele)
                bnd_chr = bnd_chr_match.group(1)

                if bnd_chr not in auto_sex_chromosomes:
                    continue
      
            bndtra_record = vcf_writer_bndtra.new_record()
            bndtra_record.chrom = record.chrom
            bndtra_record.pos = record.pos
            bndtra_record.id = f"{sv_caller_name}_{record.id}"
            bndtra_record.alleles = record.alleles
            bndtra_record.info.update(record.info)

            for i, sample in enumerate(record.samples):
                sample_data = record.samples[i]
                for key, value in sample_data.items():
                    if key == 'GQ' and value is not None:
                        bndtra_record.samples[i][key] = int(np.round(float(value)))
                    elif value is not None:
                        bndtra_record.samples[i][key] = value
      
            vcf_writer_bndtra.write(bndtra_record)
            passed_bndtra_records_written += 1
            continue
        # Skip SVs less than 50bp or have SVLEN missing
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples, changing insurveyor/survindel2 calls to dummy het genotypes except for chrY calls in the new record
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                if key == 'GT':
                    sv_record.samples[i]['GT'] = value
                    if value == (1,) and sv_record.chrom not in {"chrY", "Y"}:
                        sv_record.samples[i]['GT'] = (0, 1)
                # Round genotype qualities to nearest integer to match VCF header format
                elif key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value

        write_sv_record = True

        if 'INS' in svtype:
            sv_record.pos = sv_record.pos - 1
            sv_record.stop = sv_record.pos
            if sv_record.alts[0] == "<INS>":
                if ins_seq is not None:
                    if len(ins_seq) == svlen:
                        sv_record.alts = ((sv_record.ref + ins_seq),)
                        symbolic_ins_assembled_seq_count += 1
                    elif '-' in ins_seq:
                        ins_seq = ins_seq.replace("-", "")
                        sv_record.alts = ((sv_record.ref + ins_seq),)
                        sv_record.info['UNRESOLVED_DUPINS'] = True
                        symbolic_ins_unassembled_seq_count += 1
                    else:
                        symbolic_ins_unassembled_seq_count += 1
                else:
                    symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1

        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()
    
    ### add no svtype counter
    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_assembled_seq_count} VCF records SVTYPE is INS and symbolic. SVINSSEQ is present and fully assembled; Symbolic ALT replaced with SVINSSEQ")
    print(f"{symbolic_ins_unassembled_surveyor_seq_count} VCF records SVTYPE is INS and symbolic. SVINSSEQ is present with - separating breakend sequences. Symbolic ALT replaced with concatenated SVINSSEQ")
    print(f"{symbolic_ins_unassembled_seq_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not same length as SVLEN. Symbolic ALT kept.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records SVTYPE is INS and ALT is not symbolic")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

def process_manta_dysgu_survindel_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    symbolic_ins_unassembled_surveyor_seq_count = 0
    symbolic_ins_assembled_seq_count = 0
    symbolic_ins_unassembled_seq_count = 0
    non_symbolic_dup_count = 0
    symbolic_dup_count = 0
    symbolic_dup_assembled_seq_count = 0
    symbolic_dup_unassembled_surveyor_seq_count = 0
    symbolic_dup_unassembled_seq_count = 0
    symbolic_del_count = 0
    non_symbolic_del_count = 0
    inv_svtype_count = 0
    total_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        ins_seq = record.info.get('SVINSSEQ')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        ### add counter
        if svtype is None:
            continue
        # Skip INS with no SVLEN and Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if svlen is None:
            if 'INS' in svtype:
                continue
            if 'TRA' in svtype or 'BND' in svtype:
                if 'dysgu' in sv_caller_name:
                    continue
                alt_allele = str(record.alts[0])

                bnd_pattern = re.compile(r'[\[\]](.+?):')
                bnd_chr_match = bnd_pattern.search(alt_allele)
                bnd_chr = bnd_chr_match.group(1)

                if bnd_chr not in auto_sex_chromosomes:
                    continue
      
            bndtra_record = vcf_writer_bndtra.new_record()
            bndtra_record.chrom = record.chrom
            bndtra_record.pos = record.pos
            bndtra_record.id = f"{sv_caller_name}_{record.id}"
            bndtra_record.alleles = record.alleles
            bndtra_record.info.update(record.info)

            for i, sample in enumerate(record.samples):
                sample_data = record.samples[i]
                for key, value in sample_data.items():
                    if key == 'GQ' and value is not None:
                        bndtra_record.samples[i][key] = int(np.round(float(value)))
                    elif value is not None:
                        bndtra_record.samples[i][key] = value
      
            vcf_writer_bndtra.write(bndtra_record)
            passed_bndtra_records_written += 1
            continue
        # Skip SVs less than 50bp or have SVLEN missing
        if svlen is not None and abs(svlen) < 50:
            svlen_filter_count += 1
            continue
        # Reconstruct the SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Iterate over samples, changing insurveyor/survindel2 calls to dummy het genotypes except for chrY calls in the new record
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            for key, value in sample_data.items():
                if key == 'GT':
                    sv_record.samples[i]['GT'] = value
                    if value == (1,) and sv_record.chrom not in {"chrY", "Y"}:
                        sv_record.samples[i]['GT'] = (0, 1)
                # Round genotype qualities to nearest integer to match VCF header format
                elif key == 'GQ' and value is not None:
                    sv_record.samples[i][key] = int(np.round(float(value)))
                elif value is not None:
                    sv_record.samples[i][key] = value

        write_sv_record = True

        if 'INS' in svtype:
            sv_record.stop = sv_record.pos
            if sv_record.alts[0] == "<INS>":
                if ins_seq is not None:
                    if len(ins_seq) == svlen:
                        sv_record.alts = ((sv_record.ref + ins_seq),)
                        symbolic_ins_assembled_seq_count += 1
                    elif '-' in ins_seq:
                        ins_seq = ins_seq.replace("-", "")
                        sv_record.alts = ((sv_record.ref + ins_seq),)
                        sv_record.info['UNRESOLVED_DUPINS'] = True
                        symbolic_ins_unassembled_seq_count += 1
                    else:
                        symbolic_ins_unassembled_seq_count += 1
                else:
                    symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1

        elif 'DUP' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            sv_record.info['UNMERGED_SVTYPE'] = "DUP"
            if sv_record.alts[0] == '<DUP>' or sv_record.alts[0] == '<DUP:TANDEM>':
                if sv_record.alts[0] == '<DUP:TANDEM>':
                    sv_record.alts = ('<DUP>',)
                if ins_seq is not None:
                    if len(ins_seq) == svlen:
                        sv_record.alts = ((sv_record.ref + ins_seq),)
                        symbolic_dup_assembled_seq_count += 1
                    elif '-' in ins_seq:
                        ins_seq = ins_seq.replace("-", "")
                        sv_record.alts = ((sv_record.ref + ins_seq),)
                        sv_record.info['UNRESOLVED_DUPINS'] = True
                        symbolic_dup_unassembled_surveyor_seq_count += 1
                    else:
                        sv_record.stop = sv_record.pos + abs(svlen)
                        symbolic_dup_unassembled_seq_count += 1
                else:
                    sv_record.stop = sv_record.pos + abs(svlen)
                    symbolic_dup_count += 1
            else:
                non_symbolic_dup_count += 1

        elif 'DEL' in svtype:
            sv_record.stop = sv_record.pos + abs(svlen)
            if svlen > 0:
                sv_record.info['SVLEN'] = -svlen
            if sv_record.alts[0] == "<DEL>":
                symbolic_del_count += 1
            else:
                non_symbolic_del_count += 1
        else:
            sv_record.stop = sv_record.pos + abs(svlen)
            inv_svtype_count += 1
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()
    
    ### add no svtype counter
    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{non_auto_sex_chromosome_variants} VCF records are not chromosome 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_assembled_seq_count} VCF records SVTYPE is INS and symbolic. SVINSSEQ is present and fully assembled; Symbolic ALT replaced with SVINSSEQ")
    print(f"{symbolic_ins_unassembled_surveyor_seq_count} VCF records SVTYPE is INS and symbolic. SVINSSEQ is present with - separating breakend sequences. Symbolic ALT replaced with concatenated SVINSSEQ")
    print(f"{symbolic_ins_unassembled_seq_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not same length as SVLEN. Symbolic ALT kept.")
    print(f"{symbolic_ins_count} VCF records SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_assembled_seq_count} VCF records SVTYPE is DUP and symbolic. SVINSSEQ is present and fully assembled; Symbolic ALT replaced with SVINSSEQ")
    print(f"{symbolic_dup_unassembled_surveyor_seq_count} VCF records SVTYPE is DUP and symbolic. SVINSSEQ is present with - separating breakend sequences. Symbolic ALT replaced with concatenated SVINSSEQ")
    print(f"{symbolic_dup_unassembled_seq_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not same length as SVLEN. Symbolic ALT kept.")
    print(f"{symbolic_dup_count} VCF records SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records and SVTYPE is INV or missing.")
    print(f"{total_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Modifies VCF files output by multiple SV callers for compatibility with Truvari bench/collapse."
                    "Assumes INFO/SVLEN field value already exists in each record except for SurVeyor calls."
    )
    parser.add_argument(
        "input_vcf",
        help="Path to the input VCF file.",
    )
    parser.add_argument(
        "output_sv_vcf",
        nargs='?',
        help="Path to the output VCF file for processed SVs.",
        default=None
    )
    parser.add_argument(
        "--output-bndtra-vcf",
        help="Path to the output VCF file for BND/TRA type SVs. Required for all callers except MELT and Scramble.",
        default=None
    )
    parser.add_argument(
        "--infer-vcf-only",
        action="store_true",
        help="Infer the SV caller used to generate the VCF given as input and print the name to stdout.",
    )
    parser.add_argument(
        "--sample-name",
        default=None,
        help="New sample name to use for the output VCF.",
    )

    args = parser.parse_args()

    if not args.infer_vcf_only and args.output_sv_vcf is None:
        print(f"Error: File path for an output VCF is required.", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.input_vcf):
        print(f"Error: The input VCF file was not found: {args.input_vcf}", file=sys.stderr)
        sys.exit(1)

    # Check if VCF is MOPline headerless VCF type or Breakdancer
    mopline_headerless_vcf = False
    breakdancer_vcf = False
    try:
        with open(args.input_vcf, 'r') as f:
            lines_checked = 0
            while lines_checked < 3:
                line = f.readline()

                if lines_checked == 0 and not line.startswith('#'):
                    mopline_headerless_vcf = True
                    break
                if 'Breakdancer' in line:
                    breakdancer_vcf = True
                    break

                lines_checked += 1

    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Breakdancer VCF processing
    if breakdancer_vcf:
        if args.infer_vcf_only:
            print("breakdancer")
        
        process_breakdancer_vcf(args.input_vcf, args.output_sv_vcf, args.sample_name)
        sys.exit(0)
        
    # MOPline headerless VCF processing
    if mopline_headerless_vcf:
        if args.infer_vcf_only:
            print("mopline_headerless_vcf")
            sys.exit(0)
        
        process_mopline_headerless_vcf(args.input_vcf, args.output_sv_vcf, args.sample_name)
        sys.exit(0)

    # VCF processing for all other SV callers
    try:
        vcf_reader = pysam.VariantFile(args.input_vcf, "r")
    except ValueError as e:
        print(f"The VCF file '{args.input_vcf}' is malformed", file=sys.stderr)
        print(f"Pysam encountered an error parsing the VCF file: {e}", file=sys.stderr)
        sys.exit(1)

    first_record = None
    try:
        first_record = next(vcf_reader)
    except StopIteration:
        print(f"ERROR: VCF contains no variant records.", file=sys.stderr)
        sys.exit(1)

    vcf_reader.reset()

    if args.infer_vcf_only:
        caller = infer_input_vcf(vcf_reader.header, first_record)
        if caller:
            print(caller)
            sys.exit(0)
        else:
            print(f"ERROR: Cannot infer which SV caller generated the input VCF."
                "Valid SV caller VCFs are gcSV, INSurVeyor, Manta, SurVIndel2, dysgu, DELLY, smoove/LUMPY, MELT, Parliament2, Breakdancer, CNVnator, MOPline and GATK-SV.",
                file=sys.stderr)
            sys.exit(1)

    sv_caller_vcf = infer_input_vcf(vcf_reader.header, first_record)
    new_header = create_new_header(args.input_vcf, sv_caller_vcf, args.sample_name)

    main_sv_callers = ['manta', 'dysgu', 'survindel2']
    bndtra_output_vcf_required = main_sv_callers + ['lumpy_smoove', 'parliament', 'gcsv', 'delly', 'insurveyor']

    wham_cnvnator_callers = ['wham', 'cnvnator']

    if sv_caller_vcf in bndtra_output_vcf_required and args.output_bndtra_vcf is None:
        print(f"ERROR: An output VCF needs to be specified with --output-bndtra-vcf for the SV caller ({sv_caller_vcf}) that has been inferred from the VCF input.",
        file=sys.stderr)
        sys.exit(1)

    if sv_caller_vcf in main_sv_callers:
        process_manta_dysgu_survindel_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf in 'insurveyor':
        process_insurveyor_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf in 'gcsv':
        process_gcsv_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf in 'delly':
        process_delly_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf == 'lumpy_smoove':
        process_lumpy_smoove_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf == 'parliament':
        process_parliament_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf == 'melt':
        process_melt_vcf(vcf_reader, args.output_sv_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf in wham_cnvnator_callers:
        process_wham_cnvnator_vcf(vcf_reader, args.output_sv_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf == 'mopline':
        process_mopline_vcf(vcf_reader, args.output_sv_vcf, new_header, sv_caller_vcf)
    elif sv_caller_vcf == 'gatksv':
        process_gatksv_vcf(vcf_reader, args.output_sv_vcf, new_header, sv_caller_vcf)
    else:
        print(f"ERROR: Cannot infer which SV caller generated the input VCF."
            "Valid SV caller VCFs are gcSV, INSurVeyor, Manta, SurVIndel2, dysgu, DELLY, smoove/LUMPY, MELT, Parliament2, Breakdancer, CNVnator, MOPline and GATK-SV.",
            file=sys.stderr)
        sys.exit(1)