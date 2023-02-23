RNA-seq and ATAC-seq of influenza-specific effector memory B cells
================
Chris Scharer

# Repository Info

Scripts used for processing bulk RNA-seq and ATAC-seq data from raw fastq files through differential analysis.

**A transcriptionally distinct subset of Influenza-specific effector memory B cells predicts long-lived antibody responses to vaccination in humans**

Anoma Nellore, Esther Zumaquero, Christopher D. Scharer, Christopher F. Fucile, Christopher M. Tipton, R. Glenn King, Tian Mi, Betty Mousseau, John E. Bradley, Fen Zhou, Stuti Mutneja, Paul A. Goepfert, Jeremy M. Boss, Troy D. Randall, Ignacio Sanz, Alexander F. Rosenberg, Frances E. Lund

[Citation](https://www.biorxiv.org/content/10.1101/643973v2)

# Analysis Folders and Routines

## Tools

  - bisTools.R: common functions used for processing files


## RNA-seq

### mapping Folder

  - RNAseq.sample.manifest.txt: key file for all samples and contains
    metadata for each

  - RNAseq.pipeline.R: script for initial data organization, fastq
    qc/trimming, mapping, duplicate marking

### coverage Folder

  - geneCounts.exon.PE.R: extracts coverage for all Entrez gene
    exons and normalizes data

### differential_analysis Folder

  - diff.glm.R: pairwise differential analysis with edgeR

### Ig_coverage Folder

  - calc_Ig_coverage.R: script for calculating the coverage of each Ig segment defined in the file below from IMGT

  - IG.hg38.bed: IMGT derived bed file of Ig gene segment locations

### spliceGraph Folder

  - Igh_spliceGraph.R: script for the analysis of splicing in IgH constant regions

## ATAC-seq

### mapping Folder

  - ATACseq.sample.manifest.txt: key file for all samples and contains
    metadata for each

  - ATACseq.pipeline.bw2.R: script for initial data organization, fastq
    qc/trimming, mapping, duplicate marking, and peak calling

### coverage Folder

  - calc_coverage_.R: generates common set of peaks, extracts coverage for all peaks and normalizes data

### differential_analysis Folder

  - diff.glm.ATACseq.R: pairwise differential analysis with edgeR

