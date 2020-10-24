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


## Tutorial with example datasets
* Here is an example application with the SNU601 scDNA-seq dataset from Andor et al., 2020. 
<br/>

### Run all steps with a single command
* With the input files loaded in Step0, you can run the command for step1-6. 
```
Obj_filtered=Rundf_dna(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,
                      samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, 
                      barcodes=barcodes, size=size, assay='scDNAseq',
                      raw_counts=raw_counts, ref_counts=ref_counts, type='cellline',
                      cell_filter = 1000, SNP_filter = 20, min_vaf = 0.1, max_vaf = 0.9)
```
<br/>

#### Step0. Load the input files

* In R, load the library
```
library(Alleloscope)
```

* Set directory to the downloaded github folder.


* Read common files
```
data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
```

* Read example files for SNPs (SNP by cell matrices for ref and alt alleles)
```
barcodes=read.table("data-raw/SNU601/barcodes_sub.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readRDS("data-raw/SNU601/alt_all_sub.rds")
ref_all=readRDS("data-raw/SNU601/ref_all_sub.rds")
var_all=readRDS("data-raw/SNU601/var_all_sub.rds") 
# Info of the variants with the order the same as the order of the rows in both alt_all and ref_all
```

* Read example files for bin coverage (bin by cell matrices for tumor and normal for segmentation.)
```
raw_counts=readRDS('data-raw/SNU601/tumor_sub.rds')
ref_counts=readRDS('data-raw/SNU601/normal_sub.rds') # Normal sample from patient 6198 was used for the cell line.

# Without paired normal sample, other normal samples aligned to the same reference genome (eg. GRCh38) also work if with matched bins.
```

* Also, specify your output directory, for example:
```
dir_path <- "./data-raw/SNU601/output/"; dir.create(dir_path)
```
<br/>

#### Step1. Creating a Alleloscope object for the analysis

* First, create a Alleloscope obj
```
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scDNAseq')
```

* Filter out cells and SNPs with too few read counts
```
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=1000, SNP_filter=20, min_vaf = 0.1, max_vaf = 0.9, centro=centromere.GRCh38, telo=telomere.GRCh38) 

# suggest setting min_vaf=0.1 and max_vaf=0.9 when SNPs are called in the tumor sample for higher confident SNPs
```

* Remove obj if it takes too much space once the obj is filtered
```
#rm(Obj)
```
<br/>

#### Step2. Segmentation based on total coverage pooled across cells

* Segmentation using HMM based on total coverage

* For scATAC-seq, segmentation is suggested to be performed using matched DNA sequencing data (bulk/single).
```
Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                           raw_counts=raw_counts, # from matched DNA sequencing (bulk/single)
                           ref_counts=ref_counts, # from matched DNA sequencing (bulk/single)
                           plot_seg = TRUE)
```


![Alt text](inst/plots/segmentation.png?raw=true "SNU601 segmentation")

* Filter segments based on the numbers of SNPs.
```
Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=5000)
```

* If you want to look at specific regions, subset the seg_table_filtered
```
#(region_list=Obj_filtered$seg_table_filtered$chrr)   ## select region from this list
#sel_region_list=region_list[c(1,2)] # select region for demonstration:  "15:17760000"
#Obj_filtered$seg_table_filtered=Obj_filtered$seg_table_filtered[which(Obj_filtered$seg_table_filtered$chrr %in% sel_region_list),] # subset seg_table_filtered
## if there is no normal region in the selected regions, use segmentation plot to pick a "normal" region for normalization
```
<br/>

#### Step3. Estimate cell major haplotype proportion for each region

* Estimate theta_hat of each cell for each region in the filtered segment table (seg_table_filtered).
```
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T,cont = TRUE)

# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)
```
<br/>

#### Step4. Identify normal cells and diploid regions

* Choosing normal cells and normal regions based on the estimated theta_hat values and raw coverage across the regions.
```
Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts)

# add "select_normal" list to the Obj_filtered object. 
# The list includes barcodes for the normal cells and some candidate "normal regions"
```

* Select "normal region" from the list.
```
print(Obj_filtered$select_normal$region_normal)
Obj_filtered$ref=Obj_filtered$select_normal$region_normal[1] # choose one normal region
```
<br/>

#### Step5. Genotype each cell in each region

* Estimate cell-specific (rho_hat, theta_hat) values for each region.

* For single-cell seq data other than scDNA-seq on cell line samples, set ref_gv = "genotype_values" (from scDNA-seq) to help with rho_hat estimation.
```
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts, cov_adj=1)  # for tumor
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, ref_counts = ref_counts, cov_adj =1 ) # for cell line without normal cells in the tumor sample.
```

* Genotype all cells for each region and generate a genotype plot

* For single-cell seq data other than scDNA-seq, set ref_gt = "genotypes" (from scDNA-seq) to help with genotyping.
```
Obj_filtered=Genotype(Obj_filtered = Obj_filtered)
```
The output genotying results for two regions are shown below.

![Alt text](inst/plots/genotype.png?raw=true "SNU601 genotypes")
<br/>

#### Step6. Construct lineage structure using genotypes for each cell across all regions

* Generate lineage tree based on cell-specific genotypes across the regions.
```
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 2)
```
The output clustering result for two regions is shown below.

![Alt text](inst/plots/lineage.png?raw=true "SNU601 lineage")

* For scATAC-seq data, cells can be assigned to one of the identified subclones from matched scDNA-seq using the "AssignClones_ref" function.
<br/>

#### Save the object
```
saveRDS(Obj_filtered,paste0(dir_path, "rds/Obj_filtered.rds"))
```
<br/>

## Reference
For more information about the method, check out ##



