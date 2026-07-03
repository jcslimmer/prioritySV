# prioritySV
Merge SurVeyor, hitSV/gcSV and Manta VCFs based on their matching accuracy in specific SV size bins from benchmarking data.

## Benchmarking data and analyses
The VCFs used for benchmarking can be found at zenodo...
Scripts for preprocessing ...

## Installation
Requires Truvari v5.3.0, bcftools (only tested v1.15.1), and bedtools (only tested v2.29.0) to be available in your PATH. Also requires the pysam python library (only v0.22.1 tested). ANNOVAR must be installed and available in your PATH if using the ANNOVAR annotation option.

**Note:** Other bcftools, bedtools and pysam versions should work, but carefully check output for what is expected

Navigate to the releases page ( and download the .tar.gz archive.
