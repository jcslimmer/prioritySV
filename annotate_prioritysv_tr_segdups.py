import sys
import subprocess
import os
import tempfile
import gzip

def parse_info(info_str):
    """Parses the VCF INFO column into a dictionary."""
    info_dict = {}
    if info_str == "." or not info_str:
        return info_dict
    for item in info_str.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            info_dict[key] = val
        else:
            info_dict[item] = True
    return info_dict

def generate_query_bed(vcf_path):
    """Parses VCF and creates a temp BED file with Line Indices."""
    query_bed_fd, query_bed_path = tempfile.mkstemp(suffix=".bed")
    
    if vcf_path.endswith(".gz"):
        opener = gzip.open(vcf_path, 'rt')
    else:
        opener = open(vcf_path, 'r')

    with os.fdopen(query_bed_fd, 'w') as bed_out, opener as vcf_in:
        line_idx = 0
        for line in vcf_in:
            if line.startswith("#"):
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                line_idx += 1
                continue

            chrom = parts[0]
            try:
                pos = int(parts[1])
                ref = parts[3]
                info_col = parts[7]
                info_dict = parse_info(info_col)
                
                # Coordinate Logic
                start = pos - 1 
                if info_dict.get('SVTYPE') == 'BND':
                    end = start + 1
                elif 'END' in info_dict:
                    end = int(info_dict['END'])
                else:
                    end = start + len(ref)
                
                if start >= end: end = start + 1
                bed_out.write(f"{chrom}\t{start}\t{end}\t{line_idx}\n")
            except (ValueError, IndexError):
                pass
            line_idx += 1
            
    return query_bed_path

def calc_coverage(query_bed, ref_bed):
    """Runs bedtools coverage and returns dict {line_idx: (bp, pct)}"""
    results = {}
    cmd = ["bedtools", "coverage", "-a", query_bed, "-b", ref_bed]
    
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in process.stdout:
        cols = line.strip().split('\t')
        try:
            line_idx = int(cols[3])
            covered_bp = int(cols[5])
            fraction = float(cols[7])
            results[line_idx] = (covered_bp, fraction)
        except (ValueError, IndexError):
            continue
    process.wait()
    return results

def run_pipeline(vcf_path, bed_tr_path, bed_sd_path):
    query_bed = generate_query_bed(vcf_path)

    tr_results = calc_coverage(query_bed, bed_tr_path)

    sd_results = calc_coverage(query_bed, bed_sd_path)
    
    if vcf_path.endswith(".gz"):
        opener = gzip.open(vcf_path, 'rt')
    else:
        opener = open(vcf_path, 'r')
    
    with opener as vcf_in:
        line_idx = 0
        for line in vcf_in:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    # Add Header Definitions
                    print('##INFO=<ID=TR_OVERLAP_BP,Number=1,Type=Integer,Description="Base pairs overlapping Simple Repeats">')
                    print('##INFO=<ID=TR_OVERLAP_PCT,Number=1,Type=Float,Description="Percent overlapping Simple Repeats">')
                    print('##INFO=<ID=SD_OVERLAP_BP,Number=1,Type=Integer,Description="Base pairs overlapping Segmental Duplications">')
                    print('##INFO=<ID=SD_OVERLAP_PCT,Number=1,Type=Float,Description="Percent overlapping Segmental Duplications">')
                print(line, end='')
                continue
            
            # Retrieve Data
            tr_bp, tr_frac = tr_results.get(line_idx, (0, 0.0))
            sd_bp, sd_frac = sd_results.get(line_idx, (0, 0.0))
            
            tr_pct = round(tr_frac * 100, 2)
            sd_pct = round(sd_frac * 100, 2)

            # Construct Info String
            new_tags = []
            if tr_bp > 0:
                new_tags.append(f"TR_OVERLAP_BP={tr_bp};TR_OVERLAP_PCT={tr_pct}")
            if sd_bp > 0:
                new_tags.append(f"SD_OVERLAP_BP={sd_bp};SD_OVERLAP_PCT={sd_pct}")
            
            if new_tags:
                parts = line.strip().split('\t')
                info_col = parts[7]
                append_str = ";".join(new_tags)
                
                if info_col == "." or info_col == "":
                    parts[7] = append_str
                else:
                    parts[7] = f"{info_col};{append_str}"
                print("\t".join(parts))
            else:
                print(line, end='')
            
            line_idx += 1

    os.remove(query_bed)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python annotate_ensemblesv_tr_segdups.py <input.vcf> <simple_repeats.bed> <seg_dups.bed>")
        sys.exit(1)
    run_pipeline(sys.argv[1], sys.argv[2], sys.argv[3])
