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

SurveYor, gcSV and Manta

# Note:

**gcSV v1.0.0 was used and no longer exists due to it being renamed and re-released as HitSV (https://github.com/hitbc/HitSV). The scripts may not be compatible with potential updates made to the HitSV VCFs and I currently have no plans at the moment to update them for compatibility. gcSV v1.0.0 can be obtained from here:**
