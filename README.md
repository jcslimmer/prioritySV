# prioritySV
Merge SurVeyor, hitSV/gcSV and Manta VCFs based on benchmarking of their matching accuracy by SV type and in specific SV size bins.

For insertions/duplicatons, the priority orders in descending order per size bin are: <200bp (SurVeyor, gcSV, Manta), 200-500bp (gcSV, SurVeyor, Manta), 500bp+ (SurVeyor, gcSV, Manta).

For deletions, the priority orders in descending order per size bin are: <200bp (gcSV, SurVeyor, Manta), 200-1000bp (SurVeyor, gcSV, Manta), 1kb+ (SurVeyor, Manta, gcSV).

Inversions are merged using the priority order (SurVeyor, gcSV, Manta) if they have at least 10% reciprocal overlap and 50% size similarity. 

## Benchmarking data and analyses
The VCFs used for benchmarking (which the priority order for merging SVs is based off) can be found at (zenodo link).

Scripts for processing/standardising the VCFs for Truvari bench and running the SV database matching analyses can be found in the benchmarking_analyses folder.

## Installation
If only the SV merging is required then clone this repository:

`git clone https://github.com/jcslimmer/prioritySV`

If allele frequency estimates or tandem repeat/segmental duplication annotation of SVs are required then download the following file to the main prioritySV directory:

`wget Zenodo link with tarball`

`tar -xzvf prioritySV_afreq_trsegdup_files.tar`

## Usage

`prioritysv.sh -i gcsv.vcf,manta.vcf,surveyor.vcf -s sample_id -r reference.fa -o /path/to/output_dir`

This script processes single-sample **gcSV**, **Manta**, and **SurVeyor** VCFs. Structural variants are split by SV type and size, then merged according to the priority order defined above.

> [!NOTE]
> **Dependencies**
> - Truvari v5.3.0
> - bcftools (tested with v1.15.1)
> - bedtools (tested with v2.29.0)
> 
> **Note:** Other bcftools, bedtools and pysam versions should work, but carefully check output for any unexpected errors.
>
> All dependencies must be installed and available in your `PATH`.
>
> **Optional dependency**
> - ANNOVAR (required only when using `-a`, `--annovar-dir`)
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
>
> **Manta**
> - These scripts were developed and tested using **Manta v1.6.0**.

### Required options

| Option | Description |
| :----- | :---------- |
| `-i`, `--input-vcfs` | Comma-separated list of input SV VCF files (gcSV, Manta, SurVeyor). |
| `-s`, `--sample` | Sample name. |
| `-r`, `--ref-fasta` | Reference FASTA file. |
| `-o`, `--output-dir` | Output directory. |

### Additional options

| Option | Description |
| :----- | :---------- |
| `-h`, `--help` | Display the help message. |
| `-t`, `--threads` | Number of CPU threads to use. |
| `-a`, `--annovar-dir` | Annotate the merged SV VCF with ANNOVAR gene-based annotations. Specify the path to the `/annovar/` directory. |
| `-f`, `--allele-frequency` | Annotate the merged SV VCF with allele frequencies from gnomAD-SV, ONT-1KGP, and CoLoRS. |
| `-d`, `--repeat-sd` | Annotate the merged SV VCF with UCSC Simple Repeats and Segmental Duplications coverage. |
| `-b`, `--include-bnd` | Include BND structural variants in the merged output. |
| `-n`, `--include-nosvlen` | Include INS and INV variants with missing or zero `SVLEN` values in the merged output. |

### Allele Frequency and Repeat Annotation INFO Fields

| INFO field | Description |
| :--------- | :---------- |
| `ONT_1KGP_AF` | Allele frequency of the matched structural variant in the ONT-1KGP cohort (450 high-coverage ONT sequenced samples from the 1000 Genomes Project). |
| `ONT_1KGP_COORDS` | Coordinates of the matched ONT-1KGP structural variant. |
| `ONT_1KGP_PctSizeSim` | Size similarity (%) between the query SV and the matched ONT-1KGP SV. |
| `ONT_1KGP_PctRecOverlap` | Reciprocal overlap (%) between the query SV and the matched ONT-1KGP SV. **Note:** This is unreliable for insertions and duplications because insertions may be represented as duplications (and vice versa). It may only be useful for deletions. |
| `ONT_1KGP_SizeDiff` | Size difference (bp) between the query SV and the matched ONT-1KGP SV. |
| `ONT_1KGP_StartDist` | Difference (bp) between the start coordinates of the query SV and the matched ONT-1KGP SV. |
| `ONT_1KGP_EndDist` | Difference (bp) between the end coordinates of the query SV and the matched ONT-1KGP SV. |
| `CoLoRS_*` | Equivalent fields to `ONT_1KGP_*`, describing the matched SV in the **CoLoRS** database. |
| `gnomAD_*` | Equivalent fields to `ONT_1KGP_*`, describing the matched SV in the **gnomAD-SV** database. |
| `TR_OVERLAP_PCT` | Proportion of the SV overlapping UCSC **Simple Repeats** tracks. For example, a 10 kb deletion overlapping 2 kb of repeat sequence has a value of **20%**. Overlapping repeat intervals are merged before calculating the percentage. For insertions, this value is **100%** if the insertion breakpoint falls within a Simple Repeat. |
| `SD_OVERLAP_PCT` | Proportion of the SV overlapping UCSC **Segmental Duplications** tracks, calculated identically to `TR_OVERLAP_PCT`. |

