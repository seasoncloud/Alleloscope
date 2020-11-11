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
* GATK HaplotypeCaller (https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) is recommended to use to call germline SNPs from the standard bam files [Example script](https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/snv_calling_gatk.sh). Other SNP calling tools such as BCFtools can also be used. 
* SNPs are recommended to be called from the bam file of the matched normal samples. Without matched normal samples, our results show that calling SNPs from the tumor/cellline sample itself can also work.
<br/>
 
2. A tsv file with all cell barcodes. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/barcodes_sub.tsv)
* Each row is a barcode indicating cell identity.
* The "barcodes.tsv" files are the standard outputs of the Cell Ranger software.
<br/>
 
3. SNP by cell (sparse) matrices for both reference allele and alternative alleles. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/alt_all_sub.mtx) 
* For single-cell platforms using barcode technology with all reads in a single bam file, the VarTrix (https://github.com/10XGenomics/vartrix) tools can be used to generate SNP by cell matrices for both ref and alt alleles [Example script](https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/vartrix.sh).
* For single-cell platforms with separate bam files, the two matrices can be directly generated from multi-sample vcf files.
* The information for each SNP should be in the vcf file, the labeling for each cell should be in the barcodes.tsv file (with the same order).
<br/>
  
4. Bin by cell (sparse) matrices for tumor and normal samples. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/tumor_sub.txt) 
* The values in the matrices represent total read counts for each cell in each bin.
* Row name format:"chr1-1-20000"; The order of the columns (Each column is a cell.) should be the same as that in the barcodes.tsv.
* Without paired normal sample, other normal samples aligned to the same reference genome (eg. GRCh38) also work if with matched bins.
* CNA detection methods for scDNA-seq data such as SCOPE (https://github.com/rujinwang/SCOPE) or Cell Ranger Single Cell DNA software (https://support.10xgenomics.com/single-cell-dna/software/overview/welcome) can be used to generate the matrices for tumor and normal samples.
<br/>

## Tutorial for scDNA-seq data
* Here is an example application to the SNU601 scDNA-seq dataset from Andor et al., 2020 with five example regions. 
<br/>

### Run all steps with a single command
* With the input files loaded in Step0, you can run the command for step1-6. 
```
Obj_filtered=Rundf_dna(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,
                      samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, 
                      barcodes=barcodes, size=size, assay='scDNAseq',
                      raw_counts=raw_counts, ref_counts=ref_counts, type='cellline',
                      cell_filter = 1000, SNP_filter = 20, min_vaf = 0.1, max_vaf = 0.9,)
```
<br/>

#### Step0. Load the input files

* In R, set up the environment and read common files
```
library(Alleloscope) # load the library
setwd("~/Alleloscope/"") # set path to the github folder

dir_path <- "./samples/SNU601/scDNA/output/"; dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
```

* Read example files
```
# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw/SNU601/scDNA/barcodes_sub.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw/SNU601/scDNA/alt_all_sub.mtx")
ref_all=readMM("data-raw/SNU601/scDNA/ref_all_sub.mtx")
var_all=read.table("data-raw/SNU601/scDNA/var_all_sub.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table("data-raw/SNU601/scDNA/tumor_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)
ref_counts=read.table("data-raw/SNU601/scDNA/normal_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F) # Normal sample from patient 6198 was used for the cell line.
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

<br/>

#### Step2. Segmentation based on total coverage pooled across cells

* Segmentation using HMM based on total coverage
```
Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                           raw_counts=raw_counts, # from matched DNA sequencing (bulk/single)
                           ref_counts=ref_counts, # from matched DNA sequencing (bulk/single)
                           plot_seg = TRUE)
```


![Alt text](../../../inst/plots/segmentation.png?raw=true)
<br/>

* Filter segments based on the numbers of SNPs.
```
Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=2000)
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
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts, cov_adj=1)  # for tumor
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, ref_counts = ref_counts, cov_adj =1 ) # for cell line without normal cells in the tumor sample.
```

* Genotype all cells and generate a genotype plot for each region.
```
Obj_filtered=Genotype(Obj_filtered = Obj_filtered)
```
The output genotying results for the five regions are shown below.

![Alt text](../../../inst/plots/genotype.png?raw=true "SNU601 genotypes")
<br/><br/>

#### Step6. Construct lineage structure using genotypes for each cell across all regions

* Generate lineage tree based on cell-specific genotypes across the regions.
```
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 10)
```
The output clustering result for the five regions is shown below.

![Alt text](../../../inst/plots/lineage.png?raw=true "SNU601 lineage")
<br/><br/>

## Citation
Wu, C.-Y. et al. Alleloscope: Integrative analysis of single cell haplotype-divergent copy number alterations and chromatin accessibility changes reveals novel clonal architecture of cancers. bioRxiv (2020): [https://doi.org/10.1101/2020.10.23.349407](https://doi.org/10.1101/2020.10.23.349407)





[Back to the main page](https://github.com/seasoncloud/Alleloscope)
