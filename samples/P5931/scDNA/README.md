Alleloscope (scDNA-seq)
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description
For scDNA-seq data, Alleloscope enables allele-specific copy number profiling at the single cell level to 1. detect complex multi-allelic copy number alterations (including copy neutral loss-of-heterozygosity and mirrored subclones that have the same total copy number) and 2. reconstruct tumor lineages based on the multi-allelic copy number profile.

For more information about the method, please check out the [github](https://github.com/seasoncloud/Alleloscope) and the [paper](https://doi.org/10.1101/2020.10.23.349407).
<br/>

## Prepare input files
The following are the input files for different steps.

1. A Standard vcf file with the SNP info. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/var_all_sub.vcf)
* GATK HaplotypeCaller (https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) is recommended to use to call germline SNPs from the standard bam files ([Example script](https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/snv_calling_gatk.sh)). Other SNP calling tools such as BCFtools can also be used. 
* SNPs are recommended to be called from the bam file of the matched normal samples. Without matched normal samples, our results show that calling SNPs from the tumor/cellline sample itself can also work.
<br/>
 
2. A tsv file with all cell barcodes. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/barcodes_sub.tsv)
* Each row is a barcode indicating cell identity.
* The "barcodes.tsv" files are the standard outputs of the Cell Ranger software.
<br/>
 
3. SNP by cell (sparse) matrices for both reference allele and alternative alleles. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/alt_all_sub.mtx) 
* For single-cell platforms using barcode technology with all reads in a single bam file, the VarTrix (https://github.com/10XGenomics/vartrix) tools can be used to generate SNP by cell matrices for both ref and alt alleles ([Example script](https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/vartrix.sh)).
* For single-cell platforms with separate bam files, the two matrices can be directly generated from multi-sample vcf files.([Example script](https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/hmsns_preprocessing.sh)).
* The information for each SNP should be in the vcf file, the labeling for each cell should be in the barcodes.tsv file (with the same order).
<br/>
  
4. Bin by cell (sparse) matrices for tumor and normal samples. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/tumor_sub.txt) 
* The values in the matrices represent total read counts for each cell in each bin.
* Row name format:"chr1-1-20000"; The order of the columns (Each column is a cell.) should be the same as that in the barcodes.tsv.
* Without paired normal sample, other normal samples aligned to the same reference genome (eg. GRCh38) also work if with matched bins.
* CNA detection methods for scDNA-seq data such as SCOPE (https://github.com/rujinwang/SCOPE) or Cell Ranger Single Cell DNA software (https://support.10xgenomics.com/single-cell-dna/software/overview/welcome) can be used to generate the matrices for tumor and normal samples.
<br/>

## Tutorial for scDNA-seq data
* Here is an example application to the P5931 scDNA-seq dataset from Andor et al., 2020 with five example regions. 
<br/>

#### Step0. Load the input files

* In R, set up the environment and read common files
```
library(Alleloscope) # load the library
setwd("~/Alleloscope/") # set path to the github folder

dir_path <- "./samples/P5931/scDNA/output/"; dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
```

* Read example files
```
# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw/P5931/scDNA/barcodes_sub.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw/P5931/scDNA/alt_all_sub.mtx")
ref_all=readMM("data-raw/P5931/scDNA/ref_all_sub.mtx")
var_all=read.table("data-raw/P5931/scDNA/var_all_sub.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table("data-raw/P5931/scDNA/tumor_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)
ref_counts=read.table("data-raw/P5931/scDNA/normal_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)
```
<br/>

#### Step1. Creating a Alleloscope object for the analysis

* First, create a Alleloscope obj
```
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='P5931', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scDNAseq')
```

* Filter out cells and SNPs with too few read counts
```
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=10, SNP_filter=10, min_vaf = 0, max_vaf = 1, centro=centromere.GRCh38, telo=telomere.GRCh38) 
```

<br/>

#### Step2. Segmentation step can be skipped if using each chromosome as the segments.


#### Step3. Estimate cell major haplotype proportion for each chromosome

* Estimate theta_hat of each cell for each region in the filtered segment table (seg_table_filtered).
```
set.seed(2021)
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T, cont=T)

# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)
```
<br/>

#### Step4. Identify normal cells and diploid regions

* Choosing normal cells and normal regions based on the estimated theta_hat values and raw coverage across the regions.
```
Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts, plot_theta = TRUE)

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
```
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts)  # for tumor
```

* Genotype all cells and generate a genotype plot for each region.
```
Obj_filtered=Genotype(Obj_filtered = Obj_filtered  )
```
The output genotying results for the five regions are shown below.

![Alt text](../../../inst/plots/gtype_pre2.png?raw=true "P5931 genotypes")
<br/><br/>

#### Step6. Construct lineage structure using genotypes for each cell across all regions

* Generate lineage tree based on cell-specific genotypes across the regions.
```
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 3)
```

### (Optional) Imrpoved estimation based on the abnormal cells

* In the scatter plot of chr3 (Step5) above, there are few cells located near but not at the canonical point (1.5, 0.66). To improve the estimation, Alleloscope is able to perform second-round estimation only on these few cells with the example codes shown below. 

* Select target cells
```
sub_cells=rownames(Obj_filtered$genotypes[which(Obj_filtered$genotypes[,2]!=4),]) #select abnormal cells (4 represents diploid)
```

* 2nd-round estimation
```
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 50000, plot_stat = T, sub_cells =sub_cells, sub_region = '3', max_iter = 100 )

# theta_hat values are updated in the rds_list for the selected region.
```

* Genotype each cell in each region again
```
# Estimate genotype values
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts)

# Genotype each cell
Obj_filtered=Genotype(Obj_filtered = Obj_filtered, plot_path="./samples/P5931/scDNA/output/plots/gtype_scatter_updated.pdf")
```
The updated genotying results for chr3 are shown below.

![Alt text](../../../inst/plots/gtype_updated.png?raw=true "P5931 genotypes")
<br/><br/>



## Citation
Wu, C.-Y. et al. Alleloscope: Integrative analysis of single cell haplotype-divergent copy number alterations and chromatin accessibility changes reveals novel clonal architecture of cancers. bioRxiv (2020): [https://doi.org/10.1101/2020.10.23.349407](https://doi.org/10.1101/2020.10.23.349407)





[Back to the main page](https://github.com/seasoncloud/Alleloscope)
