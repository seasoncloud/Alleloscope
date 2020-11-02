Alleloscope 
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Tutorial for scATAC-seq data
* Here is an example application with the SU008 scATAC-seq dataset (pre-treatment) from Satpathy et al., 2019. 
<br/>

#### Step0. Load the input files

* In R, load the library
```
library(Alleloscope)
```

* Set directory to the downloaded github folder.


* Read common files
```
size=read.table("data-raw/sizes.cellranger-atac-hg19-1.2.0.txt", stringsAsFactors = F)
```

* Read example files for SNPs (SNP by cell matrices for ref and alt alleles)
```
barcodes=read.table("data-raw/SU008/scATAC/barcodes.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readRDS("data-raw/SU008/scATAC/alt_all.rds")
ref_all=readRDS("data-raw/SU008/scATAC/ref_all.rds")
var_all=readRDS("data-raw/SU008/scATAC/var_all.rds") 
# Info of the variants with the order the same as the order of the rows in both alt_all and ref_all
```

* Read example files for bin coverage (bin by cell matrices for tumor and normal for segmentation.)
```
raw_counts=readRDS('data-raw/SU008/scATAC/chr200k_fragments_sub.rds')
```

* Read known cell identity (from peaks) (optional)
```
cell_type=readRDS('data-raw/SU008/scATAC/cell_type_from_peaks.rds')
```

* Also, specify your output directory, for example:
```
dir_path <- "./samples/SU008/scATAC/output/"; dir.create(dir_path)
```
<br/>

#### Step1. Creating a Alleloscope object for the analysis

* First, create a Alleloscope obj
```
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='Sample', genome_assembly="GRCh37", dir_path=dir_path, barcodes=barcodes, size=size, assay='scATACseq')
```

* Filter out cells and SNPs with too few read counts
```
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, min_vaf = 0.1, max_vaf = 0.9) 

# suggest setting min_vaf=0.1 and max_vaf=0.9 when SNPs are called in the tumor sample for higher confident SNPs
```

* Remove obj if it takes too much space once the obj is filtered
```
#rm(Obj)
```
<br/>

#### Step2. Unbiased segmentation based on matched WES/WGS data

* Load the segmentation results of the matched WES data from Yost et al., 2019.
```
Obj_filtered$seg_table=readRDS("./data-raw/SU008/scATAC/seg_table_WES.rds")
```

* Filter segments based on the numbers of SNPs.
```
Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=500)
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

#### Step4. Retrieve a diploid region from DNA-seq data

* Assign the "normal region" from bulk DNA-seq to current scATAC-seq object
```
Obj_filtered$ref=Obj_filtered$seg_table_filtered$chrr[7] # choose one normal region
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

![Alt text](../../../inst/plots/genotype.png?raw=true "SNU601 genotypes")
<br/>

#### Step6. Construct lineage structure using cell major haplotype proportions for each cell across all regions

* Generate lineage tree based on cell-specific genotypes across the regions.
```
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 2)
```
The output clustering result for two regions is shown below.

![Alt text](../../../inst/plots/lineage.png?raw=true "SNU601 lineage")

* For scATAC-seq data, cells can be assigned to one of the identified subclones from matched scDNA-seq using the "AssignClones_ref" function.
<br/>

#### Save the object
```
saveRDS(Obj_filtered,paste0(dir_path, "rds/Obj_filtered.rds"))
```
<br/>

## Reference
For more information about the method, please check out the [github](https://github.com/seasoncloud/Alleloscope) and the [paper](https://doi.org/10.1101/2020.10.23.349407)



