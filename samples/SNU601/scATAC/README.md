Alleloscope
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Tutorial for matched scATAC-seq and scDNA-seq integration
* Here is an example application with the newly generated SNU601 scATAC-seq data and the matched scDNA-seq dataset from Andor et al., 2020.
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
barcodes=read.table("data-raw/SNU601/scATAC/barcodes.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readRDS("data-raw/SNU601/scATAC/alt_all.rds")
ref_all=readRDS("data-raw/SNU601/scATAC/ref_all.rds")
var_all=readRDS("data-raw/SNU601/scATAC/var_all.rds") 
# Info of the variants with the order the same as the order of the rows in both alt_all and ref_all
```

* Read example files for bin coverage (bin by cell matrices)
```
raw_counts=readRDS('data-raw/SNU601/scATAC/chr200k_fragments_sub.rds')
# Without paired normal sample, use matched scDNA-seq result to help normalize coverge for scATAC-seq data.
```

* Also, specify your output directory, for example:
```
dir_path <- "./samples/SNU601/scATAC/output/"; dir.create(dir_path)
```
<br/>

#### Step1. Creating a Alleloscope object for the analysis

* First, create a Alleloscope obj
```
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scATACseq')
```

* Filter out cells and SNPs with too few read counts
```
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, min_vaf = 0, max_vaf = 1, centro=centromere.GRCh38, telo=telomere.GRCh38) 
```

* Remove obj if it takes too much space once the obj is filtered
```
#rm(Obj)
```
<br/>

#### Step2. Segmentation results from matched scDNA-seq or WGS/WES

* Load information from matched scDNA-seq to assist estimation.
```
Obj_scDNA=readRDS("data-raw/SNU601/scATAC/SNU601_dna.rds")
```

* Assign segments from DNA into the current object
```
Obj_filtered$seg_table_filtered=Obj_scDNA$seg_table_filtered
```
<br/>

#### Step3. Estimate cell major haplotype proportion for each region

* Estimate theta_hat of each cell for each region in the filtered segment table (seg_table_filtered).
```
Obj_filtered = Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, min_ncell = 20, phases = Obj_scDNA$rds_list, plot_stat = T, cont = TRUE)

# The phases for each SNP estimated from DNA sequencing data can help estimate the major haplotype proportion for each cell in scATAC-seq data. 
# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)
```
<br/>

#### Step4. Retrieve a diploid region from DNA-seq data

* Assign the "normal region" from DNA-seq to current scATAC-seq object
```
Obj_filtered$ref=Obj_scDNA$ref # choose one normal region
```
<br/>

#### Step5. Genotype each cell in each region

* Estimate cell-specific (rho_hat, theta_hat) values for each region.

* Set ref_gv = "genotype_values" (from scDNA-seq) to help with rho_hat estimation for scATAC-seq data.
```
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, cov_adj =1 ,ref_gtv = Obj_scDNA$genotype_values) 
```

* Genotype all cells for each region and generate a genotype plot

* Set ref_gt = "genotypes" (from scDNA-seq) to help estimate haplotype profiles for each cell in scATAC-seq data.
```
Obj_filtered=Genotype(Obj_filtered = Obj_filtered, ref_gt = Obj_scDNA$genotypes,xmax=4)
```
The output genotying results for two regions are shown below.

![Alt text](inst/plots/genotype.png?raw=true "SNU601 genotypes")
<br/>

#### Step6. Infer clonal identity for each cell in the scATAC-seq data

* generate a clone by (marker) region matrix with values representing different haplotype profiles from matched scDNA-seq data.
```
clone.genotypes=readRDS("./data-raw/SNU601/scATAC/clone.genotypes.rds")
```

* Assign cells in scATAC-seq data to one of the subclones
```
Obj_filtered=AssignClones_ref(Obj_filtered=Obj_filtered, clone.genotypes=clone.genotypes)
```
<br/>

#### Integrate DNA-level subclones and chromatin accessibility at the single-cell level

* Perform UMAP projection using genome-wide peak profile.
```
umap_peak=readRDS("./data-raw/SNU601/scATAC/peak_umap.rds")
```

* Integrate subclones and peak signals for each cell in the scATAC-seq data
```
Clone=Obj_filtered$cloneAssign$cloneAssign[match(rownames(umap_peak), names(Obj_filtered$cloneAssign$cloneAssign))]
umap_peak=cbind(umap_peak, Clone)
```
The two signals can be visuzlized simultaneously for each cell in the scATAC-seq data. 

![Alt text](../../../inst/plots/UMAP.png?raw=true "UMAP")
<br/>

#### Save the object
```
saveRDS(Obj_filtered,paste0(dir_path, "rds/Obj_filtered.rds"))
```
<br/>

## Reference
For more information about the method, please check out the [github](https://github.com/seasoncloud/Alleloscope) and the [paper](https://doi.org/10.1101/2020.10.23.349407)



