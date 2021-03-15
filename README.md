Alleloscope 
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description
Alleloscope is a method for allele-specific copy number estimation that can be applied to single cell DNA and ATAC sequencing data (separately or in combination). Allele-specific estimation allows for the more accurate delineation of copy number states and the detection of subclonal copy-neutral loss-of-heterozygosity and mirrored CNA events. On scATAC-seq data, Alleloscope allows integrative multi-omic analysis of allele-specific copy number and chromatin accessibility for the same cell. 

For more information about the method, please check out the [paper](https://doi.org/10.1101/2020.10.23.349407).
<br/>

## Overview of Alleloscope genotyping algorithm
![Alt text](inst/plots/overview2.png?raw=true "Overview of Alleloscope genotyping algorithm")


## Install

* You can install Alleloscope with the code below:

``` R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install.packages("devtools")
devtools::install_github("seasoncloud/Alleloscope") # install
library(Alleloscope) # load
```

* You can download example datasets for Alleloscope with the following command:

  Using terminal, download the repository.
```
git clone https://github.com/seasoncloud/Alleloscope.git
```

* We have included example data in the folder data-raw/. 


## Detailed tutorials with example datasets

* Click the links below to read detailed tutorials for different data types.

1. [scDNA-seq](https://github.com/seasoncloud/Alleloscope/tree/main/samples/SNU601/scDNA)
2. [scDNA-seq with 2nd-stage estimation](https://github.com/seasoncloud/Alleloscope/tree/main/samples/P5931/scDNA)
3. [scATAC-seq](https://github.com/seasoncloud/Alleloscope/tree/main/samples/SU008/scATAC)
4. [Matched scDNA-seq and scATAC-seq](https://github.com/seasoncloud/Alleloscope/tree/main/samples/SNU601/scATAC)
<br/>

## Citation
Wu, C.-Y. et al. Alleloscope: Integrative analysis of single cell haplotype-divergent copy number alterations and chromatin accessibility changes reveals novel clonal architecture of cancers. bioRxiv (2020): [https://doi.org/10.1101/2020.10.23.349407](https://doi.org/10.1101/2020.10.23.349407)



