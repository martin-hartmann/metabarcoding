# metabarcoding

Pipelines for the analysis of metabarcoding datasets derived from Illumina paired-end sequencing. The pipelines are largely based on VSEARCH with some additional modifications. Pipelines are customized for 16S rRNA genes (commonly used for bacteria and archaea) and internal transcribed spacer regions (ITS, commonly used for fungi).

The pipelines include the following general steps.
1) Specify primers and databases
2) File formatting
3) QC
4) Quality control (phiX removal, primer trimming, PE merging, quality filtering, chimera removal, target verification)
5) Delineation into both ASVs and OTUs
6) Taxonomic classification against the UNITE database
