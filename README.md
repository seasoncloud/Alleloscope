Alleloscope 
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description
Alleloscope is a method for allele-specific copy number estimation that can be applied to single cell DNA and ATAC sequencing data (separately or in combination), allowing for integrative multi-omic analysis of allele-specific copy number and chromatin accessibility for the same cell. 

## Overview of Alleloscope genotyping algorithm
![Alt text](inst/plots/overview2.png?raw=true "Overview of Alleloscope genotyping algorithm")


## Install

* You can install Alleloscope with code below:

``` R
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
  Please download the data and store them in the same directory.


## Prepare for input files
The following are the input files for different steps.

1. A Standard vcf file with the SNP info.
* GATK HaplotypeCaller (https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) is recommended to use to call germline SNPs from the standard bam files. Other SNP calling tools such as BCFtools can also be used. 
* SNPs are recommended to be called from the bam file of the matched normal samples. Without matched normal samples, our results show that calling SNPs from the tumor/ cellline sample itself can also work. 
<br/>
 
2. A tsv file with all cell barcodes.
* Each row is a barcode indicating cell identity.
* The "barcodes.tsv" files are the standard outputs of the Cell Ranger software for different single-cell sequencing assays.
<br/>
 
3. SNP by cell matrices for both reference allele and alternative alleles. 
* For single-cell platforms using barcode technology with all reads in a single bam file, the VarTrix (https://github.com/10XGenomics/vartrix) tools can be used to generate SNP by cell matrices for both ref and alt alleles.
* For single-cell platforms with separate bam files, the two matrices can be directly generated from multi-sample vcf files.
* The information for each SNP is in the vcf file; The labeling for each cell is in the "barcodes.tsv" file (with the same order). 
<br/>
  
4. Bin by cell matrices for both tumor and normal samples. 
* The values in the matrices represent total read counts for each cell in each bin.
* Row name format:"chr1-1-20000"; Column name are different cell barcodes.
* Without paired normal sample, other normal samples aligned to the same reference genome (eg. GRCh38) also work if with matched bins.
* CNA detection methods for scDNA-seq data such as SCOPE (https://github.com/rujinwang/SCOPE) or Cell Ranger Single Cell DNA software (https://support.10xgenomics.com/single-cell-dna/software/overview/welcome) can be used to generate the matrices for tumor and normal samples.
* For 10x scATAC-seq data, peak by cell matrix can be converted to bin by cell matrix by overlaying the signals. 
<br/>

## Tutorials with example datasets
1. [scDNA-seq](https://github.com/seasoncloud/Alleloscope/tree/main/samples/SNU601/scDNA)
2. [scATAC-seq](https://github.com/seasoncloud/Alleloscope/tree/main/samples/SU008/scATAC)
3. [scDNA-seq and scATAC-seq integration](https://github.com/seasoncloud/Alleloscope/tree/main/samples/SNU601/scATAC)
<br/>

## Reference
For more information about the method, please check out the [paper](https://doi.org/10.1101/2020.10.23.349407).



