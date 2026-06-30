#!/usr/bin/env python3

import argparse
import sys
import pysam
import numpy as np
import re
import os

def infer_input_vcf(header, first_record):
    """
    Infer the SV caller VCF by searching the headers for a caller-specific identifier
    """

    header_strings = str(header)

    if 'INSurVeyor' in header_strings:
        return 'insurveyor'
    elif 'DYSGU' in header_strings:
        return 'dysgu'
    elif 'SurVIndel2' in header_strings:
        return 'survindel'
    elif 'GenerateSVCandidates' in header_strings:
        return 'manta'
    elif 'gcSV' in header_strings:
        return 'gcsv'
    elif 'LUMPY' in header_strings:
        return 'lumpy_smoove'
    elif 'MELT' in header_strings:
        return 'melt'
    elif 'WHAM' in header_strings:
        return 'wham'
    elif 'CNVnator' in header_strings:
        return 'cnvnator'
    elif 'SURVIVOR' in header_strings:
        return 'parliament'
    elif 'SurVeyor' in header_strings:
        return 'surveyor'
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

    malformed_surveyor_headers = {}
    if sv_caller == 'surveyor':
        malformed_surveyor_headers = {
            'SSP1':    ('Integer', 'Number of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1HQ':  ('Integer', 'Number of high-quality pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1mQ':  ('Integer', 'Minimum mapping quality of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1MQ':  ('Integer', 'Maximum mapping quality of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1AQ':  ('Float',   'Average mapping quality of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1SQ':  ('Float',   'Standard deviation of mapping quality of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1NMA': ('Float',   'Average NM (non-matches) of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP1NMS': ('Float',   'Standard deviation of NM (non-matches) of pairs around the breakpoint 1 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2':    ('Integer', 'Number of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2HQ':  ('Integer', 'Number of high-quality pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2mQ':  ('Integer', 'Minimum mapping quality of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2MQ':  ('Integer', 'Maximum mapping quality of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2AQ':  ('Float',   'Average mapping quality of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2SQ':  ('Float',   'Standard deviation of mapping quality of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2NMA': ('Float',   'Average NM (non-matches) of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.'),
            'SSP2NMS': ('Float',   'Standard deviation of NM (non-matches) of pairs around the breakpoint 2 in the reference allele that are discordant and yet do not support the SV.')
        }
        ids_to_replace.update(malformed_surveyor_headers.keys())
    
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
    new_header.info.add('INFERRED_SEQUENCE_SVLEN', '0', 'Flag', 'Sequence of DUP/INS was not fully assembled and no SVLEN estimate exists. SVLEN estimated from inferred sequence length.')
    new_header.info.add('AVG_SVLEN_ESTIMATE', '0', 'Flag', 'Sequence of DUP/INS was not fully assembled and or assembly not conducted. SVLEN estimate calculated from average of MINSIZE and MAXSIZE.')
    new_header.info.add('UNASSEMBLED_SVINSSEQ_SVLEN', '0', 'Flag', 'Sequence of DUP/INS was not fully assembled and or assembly not conducted. SVLEN estimate calculated from SVINSSEQ length because it was greater than average of MINSIZE and MAXSIZE.')
    new_header.info.add('LEFT_RIGHT_SVINSSEQ_SVLEN', '0', 'Flag', 'Sequence of DUP/INS was not fully assembled and no SVLEN estimate exists. SVLEN estimated from LEFT_SVINSSEQ plus RIGHT_SVINSSEQ.')

    for header_id, (header_type, header_desc) in malformed_surveyor_headers.items():
        new_header.formats.add(header_id, '.', header_type, header_desc)

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

def process_surveyor_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    pass_filter_count = 0
    non_auto_sex_chromosome_variants = 0
    hom_ref_missing_gt_count = 0
    svlen_count = 0
    svlen_filter_count = 0
    total_no_svlen_count = 0
    symbolic_no_svlen_count = 0
    non_symbolic_no_svlen_count = 0
    symbolic_assembled_seq_svlen_count = 0
    infer_ins_seq_svlen_count = 0
    no_svlen_estimated_count = 0
    N_in_svinsseq = 0
    non_symbolic_ins_count = 0
    symbolic_ins_count = 0
    non_symbolic_dup_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    non_symbolic_del_count = 0
    symbolic_inv_count = 0
    non_symbolic_inv_count = 0
    non_svtype_count = 0
    passed_bndtra_records_written = 0
    total_pass_records_written = 0

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)
    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        ins_seq = record.info.get('SVINSSEQ')
        infer_ins_seq = record.info.get('SVINSSEQ_INFERRED')
        left_ins_seq = record.info.get('LEFT_SVINSSEQ')
        right_ins_seq = record.info.get('RIGHT_SVINSSEQ')
        incomplete_assembly = record.info.get('INCOMPLETE_ASSEMBLY')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value
        first_sample = list(record.samples.keys())[0]
        gt_alleles = record.samples[first_sample].allele_indices

        # Skip SVs which are not marked as PASS
        if "PASS" not in record.samples[first_sample]['FT']:
            pass_filter_count += 1
            continue
        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        if svtype is None:
            continue
        # Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if 'TRA' in svtype or 'BND' in svtype:
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
        # Skip SVs with missing genotype for both alleles
        if all(allele is None or allele == 0 for allele in gt_alleles):
            hom_ref_missing_gt_count += 1
            continue
        # Skip SVs less than 50bp
        if svlen is not None and abs(svlen) < 50 and svlen != 0:
            svlen_filter_count += 1
            continue
        # Reconstruct SV record to associate with the new header for all SVs not filtered out that could be modified.
        sv_record = vcf_writer_sv.new_record()
        sv_record.chrom = record.chrom
        sv_record.pos = record.pos
        sv_record.id = f"{sv_caller_name}_{record.id}"
        sv_record.alleles = record.alleles
        sv_record.info.update(record.info)

        # Add ML model probability for SV existing as an integer percentage to GQ field
        for i, sample in enumerate(record.samples):
            sample_data = record.samples[i]
            epr_value = sample_data.get('EPR')
            for key, value in sample_data.items():
                if value is not None:
                    sv_record.samples[i][key] = value
            if epr_value is not None:
                epr_int_value = int(np.round(float(epr_value*100)))
                sv_record.samples[i]['GQ'] = epr_int_value
            else:
                sv_record.samples[i]['GQ'] = None

            
            # Set SVLEN to SVINSSEQ_INFERRED sequence length if it exists for insertions with no SVLEN and SVINSSEQ is unassembled
            
            final_svlen = svlen
            true_final_svlen = None
            write_sv_record = False

            # No SVLEN and ALT is symbolic 
            if final_svlen is None:
                total_no_svlen_count += 1
            if sv_record.alts[0].startswith('<') and final_svlen is None:
                symbolic_no_svlen_count += 1

                # SVINSSEQ is fully assembled
                if ins_seq is not None and not incomplete_assembly:
                    final_svlen = len(ins_seq)
                    symbolic_assembled_seq_svlen_count += 1
                # SVINSSEQ not de novo assembled, use inferred sequence for SVLEN if it exists
                elif infer_ins_seq is not None:
                    final_svlen = len(infer_ins_seq)
                    sv_record.info['INFERRED_SEQUENCE_SVLEN'] = True
                    infer_ins_seq_svlen_count += 1
                # SVINSSEQ not fully assembled and no inferred sequence exists
                else:
                    no_svlen_estimated_count += 1

            if final_svlen is not None:           
                true_final_svlen = abs(int(round(final_svlen)))

            write_sv_record = True

            if 'DEL' in svtype:
                sv_record.info['SVLEN'] = -true_final_svlen
            else:
                if true_final_svlen is not None:
                    sv_record.info['SVLEN'] = true_final_svlen

        if 'INS' in svtype:
            sv_record.stop = sv_record.pos
            if sv_record.alts[0].startswith('<'):
                symbolic_ins_count += 1
            else:
                non_symbolic_ins_count += 1
        if 'DUP' in svtype:
            sv_record.stop = sv_record.pos + abs(true_final_svlen)
            if sv_record.alts[0].startswith('<'):
                symbolic_dup_count += 1
            else:
                non_symbolic_dup_count += 1
        if 'DEL' in svtype:
            sv_record.stop = sv_record.pos + abs(true_final_svlen)
            if sv_record.alts[0].startswith('<'):
                symbolic_del_count += 1
            else:
                non_symbolic_del_count += 1
        if 'INV' in svtype:
            sv_record.stop = sv_record.pos + abs(true_final_svlen)
            if sv_record.alts[0].startswith('<'):
                symbolic_inv_count += 1
            else:
                non_symbolic_inv_count += 1
        if svtype not in ('INS', 'DUP', 'DEL', 'INV'):
            non_svtype_count += 1
            continue
        
        if write_sv_record:
            vcf_writer_sv.write(sv_record)
            total_pass_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()
    
    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{pass_filter_count} VCF records not marked as PASS. Records excluded.")
    print(f"{non_auto_sex_chromosome_variants} VCF records marked as PASS, but chromosome is not 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records marked as PASS, but SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{hom_ref_missing_gt_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but genotype is homozygous reference or missing. Records excluded.")
    print(f"{svlen_filter_count} VCF records marked as PASS, SVTYPE is not BND or TRA, genotype is non-reference, but SVLEN < 50bp. Records excluded.")
    print(f"{total_no_svlen_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but no SVLEN exists.")
    print(f"{symbolic_no_svlen_count} VCF records marked as PASS, SVTYPE is not BND or TRA, ALT is symbolic but no SVLEN exists.")
    print(f"{non_symbolic_no_svlen_count} VCF records marked as PASS, SVTYPE is not BND or TRA, ALT is not symbolic but no SVLEN exists.")
    print(f"{symbolic_assembled_seq_svlen_count} VCF records marked as PASS, SVTYPE is not BND or TRA, ALT is symbolic and no SVLEN exists, but SVINSSEQ fully assembled.")
    print(f"{infer_ins_seq_svlen_count} VCF records marked as PASS, SVTYPE is not BND or TRA, ALT is symbolic and no SVLEN exists, but SVINSSEQ is not fully assembled. Inferred sequence used for SVLEN.")
    print(f"{no_svlen_estimated_count} VCF records marked as PASS, SVTYPE is not BND or TRA, ALT is symbolic, but no SVLEN exists and it cannot be reliably estimated from any INFO fields. Records excluded.")
    print(f"{N_in_svinsseq} VCF records marked as PASS, SVTYPE is not BND or TRA, SVLEN is >= 50bp, but N's present in SVINSSEQ. Records excluded.")    
    print(f"{symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and ALT is not symbolic")
    print(f"{non_symbolic_dup_count} VCF records marked as PASS, SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is not symbolic.")
    print(f"{symbolic_inv_count} VCF records marked as PASS, SVTYPE is INV and ALT is symbolic.")
    print(f"{non_symbolic_inv_count} VCF records marked as PASS, SVTYPE is INV and ALT is not symbolic.")
    print(f"{non_svtype_count} VCF records marked as PASS and SVTYPE is not DEL/DUP/INS/INV.")
    print(f"{total_pass_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

def process_parliament_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    pass_filter_count = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    symbolic_ins_count = 0
    symbolic_dup_count = 0
    symbolic_del_count = 0
    inv_svtype_count = 0
    total_pass_records_written = 0

    new_header.formats.add('SP', '1', 'String', 'Supporting SV callers for this sample')

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    info_fields = ['SVTYPE', 'END', 'SUPP', 'CALLERS']
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = int(np.round(value[0])) if isinstance(value := record.info.get('AVGLEN'), tuple) else int(np.round(value))


        # Skip SVs which are not marked as PASS
        if 'PASS' not in record.filter and 'DUP' not in svtype:
            pass_filter_count += 1
            continue
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
            if int(record.info.get('SUPP')) < 2 or svlen >= 2000000:
                continue
            else:
                sv_record.info['SVLEN'] = svlen
                sv_record.filter.clear()
                sv_record.filter.add('PASS')
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
            total_pass_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()
    
    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{pass_filter_count} VCF records not marked as PASS. Records excluded.")
    print(f"{non_auto_sex_chromosome_variants} VCF records marked as PASS, but chromosome is not 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records marked as PASS, but SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{symbolic_dup_count} VCF records marked as PASS, SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is symbolic.")
    print(f"{inv_svtype_count} VCF records marked as PASS and SVTYPE is INV or missing.")
    print(f"{total_pass_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

def process_wham_cnvnator_vcf(vcf_reader, output_sv_vcf, new_header, sv_caller_name):

    HG38_CONTIGS = [
        ('chr1', 248956422), ('chr2', 242193529), ('chr3', 198295559),
        ('chr4', 190214555), ('chr5', 181538259), ('chr6', 170805979),
        ('chr7', 159345973), ('chr8', 145138636), ('chr9', 138394717),
        ('chr10', 133797422), ('chr11', 135086622), ('chr12', 133275309),
        ('chr13', 114364328), ('chr14', 107043718), ('chr15', 101991189),
        ('chr16', 90338345), ('chr17', 83257441), ('chr18', 80373285),
        ('chr19', 58617616), ('chr20', 64444167), ('chr21', 46709983),
        ('chr22', 50818468), ('chrX', 156040895), ('chrY', 57227415),
        ('chrM', 16569)
    ]

    for contig, length in HG38_CONTIGS:
        new_header.add_line(f"##contig=<ID={contig},length={length}>")

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    vcf_records_processed = 0
    pass_filter_count = 0
    non_auto_sex_chromosome_variants = 0
    svlen_filter_count = 0
    total_pass_records_written = 0

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    for record in vcf_reader:
        vcf_records_processed += 1
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Skip SVs which are not marked as PASS
        if 'PASS' not in record.filter:
            pass_filter_count += 1
            continue
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
        total_pass_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{pass_filter_count} VCF records not marked as PASS. Records excluded.")
    print(f"{non_auto_sex_chromosome_variants} VCF records marked as PASS, but chromosome is not 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{total_pass_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True


def process_melt_vcf(vcf_reader, output_sv_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    pass_filter_count = 0
    non_auto_sex_chromosome_variants = 0
    svlen_filter_count = 0
    total_pass_records_written = 0

    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)
    
    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Skip SVs which are not marked as PASS
        if 'PASS' not in record.filter:
            pass_filter_count += 1
            continue
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
        total_pass_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{pass_filter_count} VCF records not marked as PASS. Records excluded.")
    print(f"{non_auto_sex_chromosome_variants} VCF records marked as PASS, but chromosome is not 1-22, X or Y. Records excluded.")
    print(f"{svlen_filter_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{total_pass_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

    return True

def process_lumpy_smoove_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    pass_filter_count = 0
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
    total_pass_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Skip SVs which are not marked as PASS
        if 'PASS' not in record.filter:
            pass_filter_count += 1
            continue
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
        # Skip if N is in REF or ALT sequence
        if 'N' in record.ref or (not record.alts[0].startswith('<') and 'N' in record.alts):
            N_in_ref_alt_seq += 1
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
            total_pass_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()

    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{pass_filter_count} VCF records not marked as PASS. Records excluded.")
    print(f"{non_auto_sex_chromosome_variants} VCF records marked as PASS, but chromosome is not 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records marked as PASS, but SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_count} VCF records marked as PASS, SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records marked as PASS, SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records marked as PASS and SVTYPE is INV or missing.")
    print(f"{total_pass_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")

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

def process_manta_dysgu_insurveyor_survindel_vcf(vcf_reader, output_sv_vcf, output_bndtra_vcf, new_header, sv_caller_name):

    vcf_records_processed = 0
    pass_filter_count = 0
    non_auto_sex_chromosome_variants = 0
    passed_bndtra_records_written = 0
    svlen_filter_count = 0
    N_in_ref_alt_seq = 0
    N_in_svinsseq = 0
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
    total_pass_records_written = 0

    vcf_writer_bndtra = pysam.VariantFile(output_bndtra_vcf, "w", header=new_header)
    vcf_writer_sv = pysam.VariantFile(output_sv_vcf, "w", header=new_header)

    auto_sex_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    for record in vcf_reader:
        vcf_records_processed += 1
        svtype = record.info.get('SVTYPE')
        ins_seq = record.info.get('SVINSSEQ')
        svlen = value[0] if isinstance(value := record.info.get('SVLEN'), tuple) else value

        # Skip SVs which are not marked as PASS
        if 'PASS' not in record.filter:
            pass_filter_count += 1
            continue
        # Keep only autosome and sex chromosome calls
        if record.chrom not in auto_sex_chromosomes:
            non_auto_sex_chromosome_variants += 1
            continue
        ### add counter
        if svtype is None:
            continue
        # Skip INS with no SVLEN and Write BND/TRA type SVs to BND/TRA output VCF with new headers and GQ value for consistency between VCFs
        if svlen is None:
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
                    if 'N' in ins_seq:
                        N_in_svinsseq += 1
                        continue
                    elif len(ins_seq) == svlen:
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
                    if 'N' in ins_seq:
                        N_in_svinsseq += 1
                        continue
                    elif len(ins_seq) == svlen:
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
            total_pass_records_written += 1

    vcf_reader.close()
    vcf_writer_sv.close()
    vcf_writer_bndtra.close()
    
    ### add no svtype counter
    print(f"{vcf_records_processed} VCF records processed.")
    print(f"{pass_filter_count} VCF records not marked as PASS. Records excluded.")
    print(f"{non_auto_sex_chromosome_variants} VCF records marked as PASS, but chromosome is not 1-22, X or Y. Records excluded.")
    print(f"{passed_bndtra_records_written} VCF records marked as PASS, but SVTYPE is BND or TRA. Record written to {args.output_bndtra_vcf}")
    print(f"{svlen_filter_count} VCF records marked as PASS, SVTYPE is not BND or TRA, but SVLEN < 50bp. Records excluded.")
    print(f"{N_in_ref_alt_seq} VCF records marked as PASS, SVTYPE is not BND or TRA, SVLEN is >= 50bp, but N's present in REF or ALT sequence. Records excluded.")
    print(f"{N_in_svinsseq} VCF records marked as PASS, SVTYPE is not BND or TRA, SVLEN is >= 50bp, but N's present in SVINSSEQ. Records excluded.")    
    print(f"{symbolic_ins_assembled_seq_count} VCF records marked as PASS, SVTYPE is INS and symbolic. SVINSSEQ is present and fully assembled; Symbolic ALT replaced with SVINSSEQ")
    print(f"{symbolic_ins_unassembled_surveyor_seq_count} VCF records marked as PASS, SVTYPE is INS and symbolic. SVINSSEQ is present with - separating breakend sequences. Symbolic ALT replaced with concatenated SVINSSEQ")
    print(f"{symbolic_ins_unassembled_seq_count} VCF records marked as PASS, SVTYPE is INS and symbolic, but SVINSSEQ is not same length as SVLEN. Symbolic ALT kept.")
    print(f"{symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_ins_count} VCF records marked as PASS, SVTYPE is INS and ALT is not symbolic")
    print(f"{symbolic_dup_assembled_seq_count} VCF records marked as PASS, SVTYPE is DUP and symbolic. SVINSSEQ is present and fully assembled; Symbolic ALT replaced with SVINSSEQ")
    print(f"{symbolic_dup_unassembled_surveyor_seq_count} VCF records marked as PASS, SVTYPE is DUP and symbolic. SVINSSEQ is present with - separating breakend sequences. Symbolic ALT replaced with concatenated SVINSSEQ")
    print(f"{symbolic_dup_unassembled_seq_count} VCF records marked as PASS, SVTYPE is DUP and symbolic, but SVINSSEQ is not same length as SVLEN. Symbolic ALT kept.")
    print(f"{symbolic_dup_count} VCF records marked as PASS, SVTYPE is DUP and symbolic, but SVINSSEQ is not present. Symbolic ALT kept.")
    print(f"{non_symbolic_dup_count} VCF records marked as PASS, SVTYPE is DUP and ALT is not symbolic")
    print(f"{symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is symbolic.")
    print(f"{non_symbolic_del_count} VCF records marked as PASS, SVTYPE is DEL and ALT is not symbolic.")
    print(f"{inv_svtype_count} VCF records marked as PASS and SVTYPE is INV or missing.")
    print(f"{total_pass_records_written} VCF records written to {args.output_sv_vcf} after filters and modifications.")
    
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Modifies a VCF file using pysam: adds/verifies SVLEN INFO header, "
                    "updates SVTYPE from INFO/SIMPLE_TYPE, uses existing INFO/SVLEN for filtering, "
                    "and filters records. DUP/INS are kept if SVLEN >= 50. "
                    "Other types are kept if SVLEN <= -50 or SVLEN >= 50. "
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
            print(f"ERROR: Cannot infer which SV caller generated the input VCF. "
                "Valid SV caller VCFs are gcSV, INSurVeyor, Manta, SurVeyor, SurVIndel2, dysgu, DELLY, smoove/LUMPY, MELT, Parliament2 and CNVnator.",
                file=sys.stderr)
            sys.exit(1)

    sv_caller_vcf = infer_input_vcf(vcf_reader.header, first_record)
    new_header = create_new_header(args.input_vcf, sv_caller_vcf, args.sample_name)

    main_sv_callers = ['manta', 'dysgu', 'insurveyor', 'survindel2']
    bndtra_output_vcf_required = main_sv_callers + ['lumpy_smoove', 'parliament', 'surveyor', 'gcsv', 'delly']

    wham_cnvnator_callers = ['wham', 'cnvnator']

    if sv_caller_vcf in bndtra_output_vcf_required and args.output_bndtra_vcf is None:
        print(f"ERROR: An output VCF needs to be specified with --output-bndtra-vcf for the SV caller that has been inferred from the VCF input.",
        file=sys.stderr)
        sys.exit(1)

    if sv_caller_vcf in main_sv_callers:
        process_manta_dysgu_insurveyor_survindel_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)
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
    elif sv_caller_vcf == 'surveyor':
        process_surveyor_vcf(vcf_reader, args.output_sv_vcf, args.output_bndtra_vcf, new_header, sv_caller_vcf)

    else:
        print(f"ERROR: Cannot infer which SV caller generated the input VCF. "
            "Valid SV caller VCFs are gcSV, INSurVeyor, Manta, SurVeyor, SurVIndel2, dysgu, DELLY, smoove/LUMPY, MELT, Parliament2 and CNVnator.",
            file=sys.stderr)
        sys.exit(1)