# metabarcoding

Pipelines for the analysis of metabarcoding datasets derived from Illumina paired-end sequencing. The pipelines are largely based on VSEARCH with some additional modifications.

The pipelines include the following general steps.
1) Formatting of files to accommodate pipeline 
2) QC of raw files
3) Thorough quality control incl. phiX removal, primer trimming, PE merging, quality filtering, chimera removal, and target verification
4) Delineation into both ASVs and OTUs
5) Taxonomic classification against the UNITE database

The following script are available:
ITS2_Taylor.bash: This bash script is tailored for fungal ITS2 amplicons derived by primers 5.8S-Fun (AACTTTYRRCAAYGGATCWCT) and ITS4-Fun (AGCCTCCGCTTATTGATATGCTTAART); see https://doi.org/10.1128/AEM.02576-16 for details.
