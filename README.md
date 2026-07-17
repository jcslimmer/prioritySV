# prioritySV
Merge SurVeyor, hitSV/gcSV and Manta VCFs based on their matching accuracy in specific SV size bins from benchmarking data.

## Benchmarking data and analyses
The VCFs used for benchmarking can be found at (zenodo link)
Scripts for preprocessing ...

## Installation
Requires Truvari v5.3.0, bcftools (only tested v1.15.1), and bedtools (only tested v2.29.0) to be available in your PATH. Also requires the pysam python library (only v0.22.1 tested). ANNOVAR must be installed and available in your PATH if using the ANNOVAR annotation option.

**Note:** Other bcftools, bedtools and pysam versions should work, but carefully check output for any unexpected errors.

If only the SV merging is required then clone this repository:

`git clone https://github.com/jcslimmer/prioritySV`

If allele frequency estimates or tandem repeat/segmental duplication annotation of SVs are required then download the following file to the main prioritySV directory:

`wget Zenodo link with tarball`

`tar -xzvf prioritySV_afreq_trsegdup_files.tar`

## Usage

`prioritysv.sh -i gcsv.vcf,manta.vcf,surveyor.vcf -s sample_id -r reference.fa -o /path/to/output_dir`

This script processes single-sample **gcSV**, **Manta**, and **SurVeyor** VCFs. Structural variants are split by SV type and size, then merged according to a predefined priority order.

> [!NOTE]
> **Dependencies**
> - Truvari v5.3.0
> - bcftools (tested with v1.15.1)
> - bedtools (tested with v2.29.0)
>
> All dependencies must be installed and available in your `PATH`.
>
> **Optional dependency**
> - ANNOVAR (required only when using `--annovar-dir`)
> - The `/annovar/` directory must also be available in your `PATH`.
>
> **SV caller compatibility**
>
> **gcSV and SurVeyor**
> - These scripts were developed and tested using **gcSV v1.0.0** and **SurVeyor v0.11.0**.
> - **gcSV** has since been renamed and re-released as **HitSV**: https://github.com/hitbc/HitSV
> - The scripts may **not** be compatible with VCF format changes introduced in newer HitSV and SurVeyor releases.
> - There are currently no plans to update the scripts for compatibility with newer HitSV and SurVeyor releases. However, post an issue if running into errors and I will look into it as the fix should be fairly straightforward.
> - gcSV v1.0.0 can be obtained from **[INSERT gcSV v1.0.0 LINK HERE]**.
> - SurVeyor v0.11.0 can be obtained from **[INSERT SurVeyor v0.11.0 LINK HERE]**.

> **Manta**
> - These scripts were developed and tested using **Manta v1.6.0**.

### Required options

| Option | Description |
| :----- | :---------- |
| `-i`, `--input-vcfs` | Comma-separated list of input SV VCF files (gcSV, Manta, SurVeyor). |
| `-s`, `--sample` | Sample name. |
| `-r`, `--ref-fasta` | Reference FASTA file. |
| `-o`, `--output-dir` | Output directory. |

### Optional options

| Option | Description |
| :----- | :---------- |
| `-h`, `--help` | Display the help message. |
| `-t`, `--threads` | Number of CPU threads to use. |
| `-a`, `--annovar-dir` | Annotate the merged SV VCF with ANNOVAR gene-based annotations. Specify the path to the `/annovar/` directory. |
| `-f`, `--allele-frequency` | Annotate the merged SV VCF with allele frequencies from gnomAD-SV, ONT-1KGP, and CoLoRS. |
| `-d`, `--repeat-sd` | Annotate the merged SV VCF with UCSC Simple Repeats and Segmental Duplications coverage. |
| `-b`, `--include-bnd` | Include BND structural variants in the merged output. |
| `-n`, `--include-nosvlen` | Include INS and INV variants with missing or zero `SVLEN` values in the merged output. |
