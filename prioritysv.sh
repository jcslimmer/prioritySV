#!/usr/bin/bash -l
set -euo pipefail

### Initialise variables
script_dir="$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)"
input_vcf_array=()
allele_frequency=false
annovar_dir=""
annovar=false
include_bnd=false
include_nolen=false
repeat_sd=false
threads="1"
threads_opt=""

# Required scripts
preprocessing_script="$script_dir/truvari_merge_preprocessing.py"
af_annotation_script="$script_dir/ensemble_calls_sv_database_matching.py"
tr_sd_annotation_script="$script_dir/annotate_ensemblesv_tr_segdups.py"

# Allele frequency database VCFs
gnomad_delinv_vcf="$script_dir"/gnomad.v4.1.sv.non_neuro_controls.sites.passonly.delinv.vcf.gz
gnomad_dupins_vcf="$script_dir"/gnomad.v4.1.sv.non_neuro_controls.sites.passonly.dupins.vcf.gz
ont1kgp_delinv_vcf="$script_dir"/UWONT_450_sniffles_2.5.2_T12_merge_sorted_delinv.vcf.gz
ont1kgp_dupins_vcf="$script_dir"/UWONT_450_sniffles_2.5.2_T12_merge_sorted_dupins.vcf.gz
colors_delinv_vcf="$script_dir"/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.delinv.vcf.gz
colors_dupins_vcf="$script_dir"/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.dupins.vcf.gz

# TR and SD BED files
tr_bed="$script_dir"/simpleRepeat_filtered_noheader_noalt_named_consolidated.sorted.bed
sd_bed="$script_dir"/GRCh38_segdups_named_consolidated.sorted.bed

show_help() {
    echo "Usage: ensemble_sv_script.sh -i gcsv.vcf,manta.vcf,surveyor.vcf -s sample1 -o /path/to/output_dir"
    echo "This script processes gcSV, Manta, and SurVeyor VCFs, splits them by SV type and size then merges them based on specific priority orders."
    echo "Truvari v5.3.0, bcftools (only tested v1.15.1) and bedtools (only tested v2.29.0) must be installed and available in your PATH."
    echo "ANNOVAR must also be installed and available in your PATH (/annovar/) if using the ANNOVAR annotation option."
    echo
    echo "Required arguments:"
    echo "  -o, --output-dir    Specifies the output directory."
    echo "  -i, --input-vcfs    Specifies comma-separated VCF file paths."
    echo "  -s, --sample        Specifies the input sample name."
    echo "  -r, --ref-fasta     Specifies the reference FASTA file to use."
    echo
    echo "Optional arguments:"
    echo " -h, --help   Display help information"
    echo " -t, --threads Number of threads"
    echo " -a, --annovar-dir    Annotates merged SV VCF with ANNOVAR gene based annotations - must provide the path to the /annovar/ directory"
    echo " -f, --allele-frequency Enables gnomAD-SV, ONT-1KGP and CoLoRS allele frequency annotation of ensemble merged VCF."
    echo " -d, --repeat-sd      Annotates merged SV VCF with Simple Repeats and SegDup UCSC track coverage."
    echo " -b, --include-bnd    Include BND type SVs in the ensemble merged VCF."
    echo " -n, --include-nosvlen  Include INS/INV type SVs with missing or zero SVLEN in the ensemble merged VCF."
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
    -o|--output-dir) output_dir="$2"; shift 2 ;;
    -i|--input-vcfs) input_vcfs="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -r|--ref-fasta) ref_fasta="$2"; shift 2 ;;
    -a|--annovar-dir) annovar_dir="$2"; shift 2 ;;
    -f|--allele-frequency) allele_frequency=true; shift 1 ;;
    -d|--repeat-sd) repeat_sd=true; shift 1 ;;
    -b|--include-bnd) include_bnd=true; shift 1 ;;
    -n|--include-nosvlen) include_nosvlen=true; shift 1 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    -h|--help) show_help ;;
    *) echo "Error: Unknown option: $1"; show_help ;;
    esac
done

if [ -z "$input_vcfs" ] || [ -z "$output_dir" ] || [ -z "$sample" ] || [ -z "$ref_fasta" ]; then
    echo "Error: --output-dir --input-vcfs --sample and --ref-fasta must be provided."
    show_help
fi

if [[ "$threads" =~ ^[0-9]+$ ]]; then
    threads_opt="--threads $threads"
else
    echo "Error: The value for --threads must be a positive integer."
    show_help
fi

IFS="," read -r -a input_vcf_array <<< "$input_vcfs"

if [ ${#input_vcf_array[@]} -ne 3 ]; then
    echo "Error: Three comma-separated input VCFs are required: gcSV, Manta, SurVeyor"
    show_help
fi

### Check required ref fasta, script and tools
if [ ! -f "$ref_fasta" ]; then
    echo "Error: The reference FASTA file was not found: $ref_fasta"
    exit 1
fi

if [ ! -f "$preprocessing_script" ]; then
    echo "Error: The required Python script 'truvari_merge_preprocessing.py' was not found in $script_dir"
    exit 1
fi

for tool in bcftools truvari bedtools bgzip tabix; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool not found in PATH."
        exit 1
    fi
done

# Validate optional allele frequency and gene annotation argument files/paths
if [ "$allele_frequency" = true ]; then
    if [ ! -f "$af_annotation_script" ]; then
        echo "Error: The required Python script 'ensemble_calls_sv_database_matching.py' was not found in $script_dir"
        exit 1
    fi
    for db_vcf in "$ont1kgp_delinv_vcf" "$ont1kgp_dupins_vcf" "$colors_delinv_vcf" "$colors_dupins_vcf" "$gnomad_delinv_vcf" "$gnomad_dupins_vcf"; do
        if [ ! -f "$db_vcf" ]; then
            echo "Error: Allele frequency annotation enabled but required VCFs were not found: $db_vcf"
            exit 1
        fi
    done
fi

if [ "$repeat_sd" = true ]; then
    if [ ! -f "$tr_sd_annotation_script" ]; then
        echo "Error: The required Python script 'annotate_ensemblesv_tr_segdups.py' was not found in $script_dir"
        exit 1
    fi
    if [ ! -f "$tr_bed" ] || [ ! -f "$sd_bed" ]; then
        echo "Error: Required BED files for Simple Repeat and Segmental Duplication annotation were not found."
        exit 1
    fi
fi

if [ -n "$annovar_dir" ]; then
    if [ ! -d "$annovar_dir" ]; then
        echo "Error: The ANNOVAR directory provided via --annovar-dir was not found:"
        echo "$annovar_dir"
        exit 1
    fi

    table_annovar_script="$annovar_dir/table_annovar.pl"
    humandb_dir="$annovar_dir/humandb"

    if [ ! -x "$table_annovar_script" ]; then
        echo "Error: The 'table_annovar.pl' script was not found or is not executable in the ANNOVAR directory."
        echo "$table_annovar_script"
        exit 1
    fi

    if [ ! -d "$humandb_dir" ]; then
        echo "Error: The 'humandb' directory was not found in the ANNOVAR directory."
        echo "$humandb_dir"
        exit 1
    fi

    ensgene_txt_file="$humandb_dir/hg38_ensGene.txt"
    ensgene_fa_file="$humandb_dir/hg38_ensGeneMrna.fa"

    if [ ! -f "$ensgene_txt_file" ] || [ ! -f "$ensgene_fa_file" ]; then
        echo "Error: The ensGene database files for hg38 were not found in the 'humandb' directory."
        echo "$ensgene_txt_file"
        echo "$ensgene_fa_file"
        exit 1
    fi
    annovar=true
fi

declare -A sv_caller_vcf

for vcf_file in "${input_vcf_array[@]}"; do
    if [ ! -f "$vcf_file" ]; then
        echo "Error: Input VCF not found: '$vcf_file'"
        exit 1
    fi
    caller_name=$(python3 "$preprocessing_script" "$vcf_file" --infer-vcf-only)

    if [ -z "$caller_name" ] || [ "$caller_name" == *"ERROR"* ]; then
        echo "Error: Could not infer SV caller VCF from '$vcf_file'."
        exit 1
    fi
    sv_caller_vcf["$caller_name"]="$vcf_file"
done

if [ -z "${sv_caller_vcf[gcsv]:-}" ] ||[ -z "${sv_caller_vcf[manta]:-}" ] ||[ -z "${sv_caller_vcf[surveyor]:-}" ]; then
    echo "Error: Could not match inputs to required callers (gcSV, Manta, SurVeyor)."
    echo "Found: ${!sv_caller_vcf[@]}"
    exit 1
fi

mkdir -p "$output_dir"/workdir/{gcsv,manta,surveyor}

gcsv_output_dir="$output_dir"/workdir/gcsv
manta_output_dir="$output_dir"/workdir/manta
surveyor_output_dir="$output_dir"/workdir/surveyor

gcsv_processed_prefix="$sample"_gcsv_truvari_preprocessing
manta_processed_prefix="$sample"_manta_truvari_preprocessing
surveyor_processed_prefix="$sample"_surveyor_truvari_preprocessing

mkdir -p "$output_dir"/workdir/truvari_merge/

output_truvari_dir="$output_dir"/workdir/truvari_merge/

### preprocess SV caller VCFs
python3 "$preprocessing_script" \
"${sv_caller_vcf[gcsv]}" \
"$gcsv_output_dir"/"$gcsv_processed_prefix".vcf.gz \
--output-bndtra-vcf "$gcsv_output_dir"/"$gcsv_processed_prefix"_bndtra.vcf.gz \
--sample-name "$sample"

python3 "$preprocessing_script" \
"${sv_caller_vcf[surveyor]}" \
"$surveyor_output_dir"/"$surveyor_processed_prefix".vcf.gz \
--output-bndtra-vcf "$surveyor_output_dir"/"$surveyor_processed_prefix"_bndtra.vcf.gz \
--sample-name "$sample"

python3 "$preprocessing_script" \
"${sv_caller_vcf[manta]}" \
"$manta_output_dir"/"$manta_processed_prefix".vcf.gz \
--output-bndtra-vcf "$manta_output_dir"/"$manta_processed_prefix"_bndtra.vcf.gz \
--sample-name "$sample"

if [ -n "$threads" ]; then
    threads_opt="--threads $threads"
fi

### extract INS/INV SVs with no or zero SVLEN and/or BND VCFs if enabled
nosvlen_vcf_list=()
bndtra_vcf_list=()

for prefix in "$gcsv_output_dir/$gcsv_processed_prefix" "$manta_output_dir/$manta_processed_prefix" "$surveyor_output_dir/$surveyor_processed_prefix"; do
    
    bcftools sort "${prefix}.vcf.gz" -Oz -o "${prefix}.sorted.vcf.gz"
    mv "${prefix}.sorted.vcf.gz" "${prefix}.vcf.gz"

    if [ "$include_nosvlen" = true ]; then
        bcftools view -i '(INFO/SVTYPE="INS" && (INFO/SVLEN="." || INFO/SVLEN==0)) || (INFO/SVTYPE="INV" && INFO/SVLEN==0)' \
            $threads_opt "${prefix}.vcf.gz" -Oz -o "${prefix}_nosvlen.vcf.gz"
        
        tabix -p vcf "${prefix}_nosvlen.vcf.gz"
        
        nosvlen_vcf_list+=("${prefix}_nosvlen.vcf.gz")
    fi

    if [ "$include_bnd" = true ]; then
        bcftools sort "${prefix}_bndtra.vcf.gz" -Oz -o "${prefix}_bndtra.sorted.vcf.gz"
        tabix -p vcf "${prefix}_bndtra.sorted.vcf.gz"
        
        bndtra_vcf_list+=("${prefix}_bndtra.sorted.vcf.gz")
    fi
done

### concatenate all callers nosvlen and/or BND VCFs if enabled
if [ "$include_nosvlen" = true ]; then
    all_nosvlen_vcf="$output_truvari_dir/${sample}_all_callers_nosvlen.vcf.gz"
    bcftools concat $threads_opt --allow-overlaps "${nosvlen_vcf_list[@]}" -Oz -o "$all_nosvlen_vcf"
    bcftools sort "$all_nosvlen_vcf" -Oz -o "${all_nosvlen_vcf%.vcf.gz}.sorted.vcf.gz"
    tabix -p vcf "${all_nosvlen_vcf%.vcf.gz}.sorted.vcf.gz"
    all_nosvlen_vcf="${all_nosvlen_vcf%.vcf.gz}.sorted.vcf.gz"
fi

if [ "$include_bnd" = true ]; then
    all_bndtra_vcf="$output_truvari_dir/${sample}_all_callers_bndtra.vcf.gz"
    bcftools concat $threads_opt --allow-overlaps "${bndtra_vcf_list[@]}" -Oz -o "$all_bndtra_vcf"
    bcftools sort "$all_bndtra_vcf" -Oz -o "${all_bndtra_vcf%.vcf.gz}.sorted.vcf.gz"
    tabix -p vcf "${all_bndtra_vcf%.vcf.gz}.sorted.vcf.gz"
    all_bndtra_vcf="${all_bndtra_vcf%.vcf.gz}.sorted.vcf.gz"
fi

### function for splitting SV caller VCFs by size and type
split_caller_vcf() {
    local input_vcf="$1"
    local output_base="$2"
    
    # small INS/DUP < 200bp
    bcftools view -i '(INFO/SVTYPE="DUP" || INFO/SVTYPE="INS") && ABS(INFO/SVLEN) < 200' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_dupins_small.vcf.gz"
    
    # medium INS/DUP 200-500bp
    bcftools view -i '(INFO/SVTYPE="DUP" || INFO/SVTYPE="INS") && ABS(INFO/SVLEN) >= 200 && ABS(INFO/SVLEN) <= 500' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_dupins_medium.vcf.gz"
    
    # large INS/DUP > 500bp
    bcftools view -i '(INFO/SVTYPE="DUP" || INFO/SVTYPE="INS") && ABS(INFO/SVLEN) > 500' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_dupins_large.vcf.gz"
    
    # small DEL < 200bp
    bcftools view -i 'INFO/SVTYPE="DEL" && ABS(INFO/SVLEN) < 200' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_del_small.vcf.gz"
    
    # medium DEL 200-1000bp
    bcftools view -i 'INFO/SVTYPE="DEL" && ABS(INFO/SVLEN) >= 200 && ABS(INFO/SVLEN) <= 1000' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_del_medium.vcf.gz"

    # large DEL > 1000bp
    bcftools view -i 'INFO/SVTYPE="DEL" && ABS(INFO/SVLEN) > 1000' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_del_large.vcf.gz"

    # INV
    bcftools view -i 'INFO/SVTYPE="INV"' \
        $threads_opt "$input_vcf" -Oz -o "${output_base}_inv.vcf.gz"
    
    # index split VCFs
    for f in "${output_base}"_*.vcf.gz; do 
        if [[ "$f" == *"bndtra"* ]] || [[ "$f" == *"nosvlen"* ]]; then continue; fi
        tabix -p vcf "$f"
    done
}

split_caller_vcf "$gcsv_output_dir/$gcsv_processed_prefix.vcf.gz" "$gcsv_output_dir/$gcsv_processed_prefix"
split_caller_vcf "$manta_output_dir/$manta_processed_prefix.vcf.gz" "$manta_output_dir/$manta_processed_prefix"
split_caller_vcf "$surveyor_output_dir/$surveyor_processed_prefix.vcf.gz" "$surveyor_output_dir/$surveyor_processed_prefix"

### function for merging using truvari bench
# runs truvari bench, combines the false positive calls from fp.vcf with tp-base.vcf to create a new base VCF for further merging with bench

run_priority_merge() {
    local svtype_size="$1"
    local truvari_params="$2"
    shift 2
    local ordered_vcfs=("$@") # array of VCFs in priority order (determined by the order they are passed to the function)
     
    # initialise highest priority caller VCF
    local first_base_vcf="${ordered_vcfs[0]}"
    local current_base_vcf="$output_truvari_dir/${svtype_size}_merged_step0.vcf.gz"
    
    # copy highest priority caller VCF to temp VCF as initial base VCF
    cp "$first_base_vcf" "$current_base_vcf"
    tabix -p vcf "$current_base_vcf"
    
    # iterate through caller VCFs in priority order
    for (( i=1; i<${#ordered_vcfs[@]}; i++ )); do
        local comp_vcf="${ordered_vcfs[$i]}"
        local priority_step_dir="$output_truvari_dir/${svtype_size}_step${i}_bench"
        
        truvari bench \
            -f "$ref_fasta" \
            -b "$current_base_vcf" \
            -c "$comp_vcf" \
            $truvari_params \
            -o "$priority_step_dir"
        
        # new priority merged step VCF = current base VCF + unique SV calls from current priority merging step (fp.vcf.gz)
        local priority_merged_step_vcf="$output_truvari_dir/${svtype_size}_merged_step${i}.vcf.gz"
        
        bcftools concat \
            $threads_opt \
            --allow-overlaps \
            "$current_base_vcf" \
            "$priority_step_dir/fp.vcf.gz" \
            -Oz -o "$priority_merged_step_vcf"
            
        tabix -p vcf "$priority_merged_step_vcf"
        current_base_vcf="$priority_merged_step_vcf"
    done
    
    # move and rename final merged output VCF for current SV type and size bin
    mv "$current_base_vcf" "$output_truvari_dir/${sample}_${svtype_size}_final_merged.vcf.gz"
    tabix -p vcf "$output_truvari_dir/${sample}_${svtype_size}_final_merged.vcf.gz"
}

# small INS/DUP < 200bp
# Priority: SurVeyor -> gcSV -> Manta
run_priority_merge "dupins_small" \
    "--pick single --sizemax 10000000 --pctsize 0.7 --refdist 500 --dup-to-ins -p 0" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_dupins_small.vcf.gz" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_dupins_small.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_dupins_small.vcf.gz"

# medium INS/DUP 200-500bp
# Priority: gcSV -> SurVeyor -> Manta
run_priority_merge "dupins_medium" \
    "--pick single --sizemax 10000000 --pctsize 0.7 --refdist 500 --dup-to-ins -p 0" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_dupins_medium.vcf.gz" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_dupins_medium.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_dupins_medium.vcf.gz"

# large INS/DUP > 500bp
# Priority: SurVeyor -> gcSV -> Manta
run_priority_merge "dupins_large" \
    "--pick single --sizemax 10000000 --pctsize 0.7 --refdist 500 --dup-to-ins -p 0" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_dupins_large.vcf.gz" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_dupins_large.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_dupins_large.vcf.gz"

# small DEL < 200bp
# Priority: gcSV -> SurVeyor -> Manta
run_priority_merge "del_small" \
    "--pick single --sizemax 10000000 --pctovl 0.7 --pctsize 0.7 -p 0" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_del_small.vcf.gz" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_del_small.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_del_small.vcf.gz"

# medium DEL 200-1000bp
# Priority: SurVeyor -> gcSV -> Manta
run_priority_merge "del_medium" \
    "--pick single --sizemax 10000000 --pctovl 0.7 --pctsize 0.7 -p 0" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_del_medium.vcf.gz" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_del_medium.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_del_medium.vcf.gz"

# large DEL > 1kb
# Priority: SurVeyor -> Manta -> gcSV
run_priority_merge "del_large" \
    "--pick single --sizemax 10000000 --pctovl 0.7 --pctsize 0.7 -p 0" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_del_large.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_del_large.vcf.gz" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_del_large.vcf.gz"

# Inversions
run_priority_merge "inv" \
    "--pick single --sizemax 10000000 --pctovl 0.1 --pctsize 0.5 -p 0" \
    "$surveyor_output_dir/${surveyor_processed_prefix}_inv.vcf.gz" \
    "$gcsv_output_dir/${gcsv_processed_prefix}_inv.vcf.gz" \
    "$manta_output_dir/${manta_processed_prefix}_inv.vcf.gz"

# concatenate DUP/INS priority merged VCFs
bcftools concat \
    $threads_opt \
    --allow-overlaps \
    "$output_truvari_dir/${sample}_dupins_small_final_merged.vcf.gz" \
    "$output_truvari_dir/${sample}_dupins_medium_final_merged.vcf.gz" \
    "$output_truvari_dir/${sample}_dupins_large_final_merged.vcf.gz" \
    -Oz -o "$output_truvari_dir/${sample}_dupins_priority_merged.vcf.gz"

bcftools sort \
    "$output_truvari_dir/${sample}_dupins_priority_merged.vcf.gz" \
    -Oz -o "$output_truvari_dir/${sample}_dupins_priority_merged.sorted.vcf.gz"
tabix -p vcf "$output_truvari_dir/${sample}_dupins_priority_merged.sorted.vcf.gz"

# oncatenate DEL and INV priority merged VCFs
bcftools concat \
    $threads_opt \
    --allow-overlaps \
    "$output_truvari_dir/${sample}_del_small_final_merged.vcf.gz" \
    "$output_truvari_dir/${sample}_del_medium_final_merged.vcf.gz" \
    "$output_truvari_dir/${sample}_del_large_final_merged.vcf.gz" \
    "$output_truvari_dir/${sample}_inv_final_merged.vcf.gz" \
    -Oz -o "$output_truvari_dir/${sample}_delinv_priority_merged.vcf.gz"

bcftools sort \
    "$output_truvari_dir/${sample}_delinv_priority_merged.vcf.gz" \
    -Oz -o "$output_truvari_dir/${sample}_delinv_priority_merged.sorted.vcf.gz"
tabix -p vcf "$output_truvari_dir/${sample}_delinv_priority_merged.sorted.vcf.gz"

# current VCF variables for allele frequency and ANNOVAR options
current_dupins_vcf="$output_truvari_dir/${sample}_dupins_priority_merged.sorted.vcf.gz"
current_delinv_vcf="$output_truvari_dir/${sample}_delinv_priority_merged.sorted.vcf.gz"

### allele frequency annotation
if [ "$allele_frequency" = true ]; then

    # function to run allele freq annotation
    # args: <database_vcf> <input_vcf> <output_dir_suffix> <output_vcf_name> <svtype_truvari_args>
    run_af_annotation() {
        local database_vcf=$1
        local input_vcf=$2
        local out_dir=$3
        local out_vcf=$4
        local args=$5
      
        truvari bench \
            -b "$database_vcf" \
            -f "$ref_fasta" \
            -c "$input_vcf" \
            --pick single \
            -S 50 \
            -s 50 \
            --passonly \
            --sizemax 10000000 \
            $args \
            -p 0 \
            -o "$out_dir"

        python3 "$af_annotation_script" \
            -b "$out_dir/tp-base.vcf.gz" \
            -c "$out_dir/tp-comp.vcf.gz" \
            -f "$out_dir/fp.vcf.gz" \
            -o "$out_vcf"

        bcftools sort "$out_vcf" -Oz -o "${out_vcf%.vcf.gz}.sorted.vcf.gz"
        tabix -p vcf "${out_vcf%.vcf.gz}.sorted.vcf.gz"
    }

    # ONT-1KGP
    run_af_annotation "$ont1kgp_delinv_vcf" "$current_delinv_vcf" \
        "$output_truvari_dir/ONT1KGP_delinv_pctovl70" \
        "$output_truvari_dir/${sample}_delinv_ont1kgp.vcf.gz" \
        "--pctovl 0.7"

    run_af_annotation "$ont1kgp_dupins_vcf" "$current_dupins_vcf" \
        "$output_truvari_dir/ONT1KGP_dupins_refdist500_pctsize70" \
        "$output_truvari_dir/${sample}_dupins_ont1kgp.vcf.gz" \
        "--refdist 500 --dup-to-ins --pctsize 0.7"

    # update current VCFs
    current_delinv_vcf="$output_truvari_dir/${sample}_delinv_ont1kgp.sorted.vcf.gz"
    current_dupins_vcf="$output_truvari_dir/${sample}_dupins_ont1kgp.sorted.vcf.gz"

    # CoLoRS
    run_af_annotation "$colors_delinv_vcf" "$current_delinv_vcf" \
        "$output_truvari_dir/CoLoRS_delinv_pctovl70" \
        "$output_truvari_dir/${sample}_delinv_ont1kgp_colors.vcf.gz" \
        "--pctovl 0.7"

    run_af_annotation "$colors_dupins_vcf" "$current_dupins_vcf" \
        "$output_truvari_dir/CoLoRS_dupins_refdist500_pctsize70" \
        "$output_truvari_dir/${sample}_dupins_ont1kgp_colors.vcf.gz" \
        "--refdist 500 --dup-to-ins --pctsize 0.7"

    # update current VCFs
    current_delinv_vcf="$output_truvari_dir/${sample}_delinv_ont1kgp_colors.sorted.vcf.gz"
    current_dupins_vcf="$output_truvari_dir/${sample}_dupins_ont1kgp_colors.sorted.vcf.gz"

    # gnomAD
    run_af_annotation "$gnomad_delinv_vcf" "$current_delinv_vcf" \
        "$output_truvari_dir/gnomAD_delinv_pctovl70" \
        "$output_truvari_dir/${sample}_delinv_ont1kgp_colors_gnomad.vcf.gz" \
        "--pctovl 0.7"

    run_af_annotation "$gnomad_dupins_vcf" "$current_dupins_vcf" \
        "$output_truvari_dir/gnomAD_dupins_refdist500_pctsize70" \
        "$output_truvari_dir/${sample}_dupins_ont1kgp_colors_gnomad.vcf.gz" \
        "--refdist 500 --dup-to-ins --pctsize 0.7"

    # update current VCFs
    current_delinv_vcf="$output_truvari_dir/${sample}_delinv_ont1kgp_colors_gnomad.sorted.vcf.gz"
    current_dupins_vcf="$output_truvari_dir/${sample}_dupins_ont1kgp_colors_gnomad.sorted.vcf.gz"
fi

### Merge final delinv and dupins VCFs
merged_dupins_delinv_vcf="$output_truvari_dir/${sample}_surveyor_manta_gcsv_ensemble_merged.vcf.gz"

bcftools concat \
    $threads_opt \
    --allow-overlaps \
    "$current_dupins_vcf" \
    "$current_delinv_vcf" \
    -Oz -o "$merged_dupins_delinv_vcf"

bcftools sort "$merged_dupins_delinv_vcf" -Oz -o "${merged_dupins_delinv_vcf%.vcf.gz}.sorted.vcf.gz"
tabix -p vcf "${merged_dupins_delinv_vcf%.vcf.gz}.sorted.vcf.gz"

#main VCF pointer
current_main_vcf="${merged_dupins_delinv_vcf%.vcf.gz}.sorted.vcf.gz"

### Filter out BNDs and INS/INVs with no SVLEN overlapping priority merged SVs (10bp) then add them to final VCF
# only done if -b and/or -n options used
if [ "$include_nosvlen" = true ]; then
    filtered_nosvlen_vcf="$output_truvari_dir/${sample}_nosvlen_overlap_filtered.vcf.gz"
    bedtools window -v -w 10 -header \
    -a <(bcftools view $threads_opt "$all_nosvlen_vcf") \
    -b <(bcftools view $threads_opt "$current_main_vcf") | bgzip > "$filtered_nosvlen_vcf"
    tabix -p vcf "$filtered_nosvlen_vcf"

    nosvlen_appended_vcf="${current_main_vcf%.sorted.vcf.gz}_with_nosvlen.vcf.gz"
    bcftools concat $threads_opt --allow-overlaps \
    "$current_main_vcf" \
    "$filtered_nosvlen_vcf" \
    -Oz -o "$nosvlen_appended_vcf"
    bcftools sort "$nosvlen_appended_vcf" -Oz -o "${nosvlen_appended_vcf%.vcf.gz}.sorted.vcf.gz"
    tabix -p vcf "${nosvlen_appended_vcf%.vcf.gz}.sorted.vcf.gz"
    #update main VCF pointer
    current_main_vcf="${nosvlen_appended_vcf%.vcf.gz}.sorted.vcf.gz"
fi

if [ "$include_bnd" = true ]; then
    filtered_bnd_vcf="$output_truvari_dir/${sample}_bndtra_overlap_filtered.vcf.gz"
    bedtools window -v -w 10 -header \
    -a <(bcftools view $threads_opt "$all_bndtra_vcf") \
    -b <(bcftools view $threads_opt "$current_main_vcf") | bgzip > "$filtered_bnd_vcf"
    tabix -p vcf "$filtered_bnd_vcf"

    bnd_appended_vcf="${current_main_vcf%.sorted.vcf.gz}_with_bndtra.vcf.gz"
    
    bcftools concat $threads_opt --allow-overlaps \
    "$current_main_vcf" \
    "$filtered_bnd_vcf" \
    -Oz -o "$bnd_appended_vcf"
    bcftools sort "$bnd_appended_vcf" -Oz -o "${bnd_appended_vcf%.vcf.gz}.sorted.vcf.gz"
    tabix -p vcf "${bnd_appended_vcf%.vcf.gz}.sorted.vcf.gz"
    #update main VCF pointer
    current_main_vcf="${bnd_appended_vcf%.vcf.gz}.sorted.vcf.gz"
fi

### repeat and segdup annotation
if [ "$repeat_sd" = true ]; then
    tr_sd_annotated_vcf="$output_truvari_dir/${sample}_tr_sd_annotated.vcf.gz"
    
    python3 "$tr_sd_annotation_script" "$current_main_vcf" "$tr_bed" "$sd_bed" | bgzip -c -@ "$threads" > "$tr_sd_annotated_vcf"
    tabix -p vcf "$tr_sd_annotated_vcf"
    
    # update main VCF pointer
    current_main_vcf="$tr_sd_annotated_vcf"
fi

### ANNOVAR anotation
if [ "$annovar" = true ]; then
    
    mkdir -p "$output_dir"/workdir/annovar/

    input_basename=$(basename "$current_main_vcf" .vcf.gz)
    annovar_output_prefix="$output_dir/workdir/annovar/${input_basename}_annovar"
    
    "$annovar_dir"/table_annovar.pl \
        "$current_main_vcf" \
        "$annovar_dir/humandb/" \
        -buildver hg38 \
        -out "$annovar_output_prefix" \
        -protocol ensGene \
        -operation g \
        -vcfinput \
        -remove \
        -thread "$threads" \
        > "$output_dir/workdir/annovar/${sample}_annovar_annotation.log" 2>&1

    annovar_vcf="${annovar_output_prefix}.hg38_multianno.vcf"
    
    bgzip -c "$annovar_vcf" > "${annovar_vcf}.gz"
    tabix -p vcf "${annovar_vcf}.gz"
    # update main VCF pointer
    current_main_vcf="${annovar_vcf}.gz"
fi

# final VCF sorting and naming
final_vcf_name="${sample}_surveyor_gcsv_manta_ensemble_merged"
annotation_suffix=""

if [ "$allele_frequency" = true ]; then annotation_suffix="${annotation_suffix}_allele_frequency"; fi
if [ "$annovar" = true ]; then annotation_suffix="${annotation_suffix}_annovar"; fi
if [ "$repeat_sd" = true ]; then annotation_suffix="${annotation_suffix}_repeat"; fi

if [ -n "$annotation_suffix" ]; then
    final_output_vcf="${final_vcf_name}${annotation_suffix}_annotated.sorted.vcf.gz"
else
    final_output_vcf="${final_vcf_name}.vcf.gz"
fi

bcftools sort "$current_main_vcf" -Oz -o "$output_dir/${final_output_vcf}"
tabix -p vcf "$output_dir/${final_output_vcf}"