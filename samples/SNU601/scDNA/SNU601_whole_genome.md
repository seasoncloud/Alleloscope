Alleloscope
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania


#### Step0. Load the input files

* In R, set up the environment and read common files
```
library(Alleloscope) # load the library

dir_path <- "~/alleloscope_output/"; dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("~/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
```

* Read example files
```
# SNP by cell matrices for ref and alt alleles
barcodes=readRDS("~/barcodes.rds)
alt_all=readRDS("~/alt_all.rds")
ref_all=readRDS("~/ref_all.rds")
var_all=readRDS("~/var_all.rds")

# bin by cell matrices for tumor and normal for segmentation
raw_counts=readRDS("~/raw_counts.rds")
ref_counts=readRDS("~/ref_counts_from_P6198.rds")
Normal sample from patient 6198 was used for the cell line.
```
<br/>

#### Step1. Creating a Alleloscope object for the analysis

* First, create a Alleloscope obj
```
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='SNU601', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scDNAseq')
```

* Filter out cells and SNPs with too few read counts
```
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=1000, SNP_filter=20, min_vaf = 0, max_vaf = 1, centro=centromere.GRCh38, telo=telomere.GRCh38) 

```

<br/>

#### Step2. Segmentation based on total coverage pooled across cells

* Segmentation using HMM based on total coverage
```
Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                           raw_counts=raw_counts, # from scDNA-seq sequencing 
                           ref_counts=ref_counts, # from other scDNA-seq sequencing
                           hmm_states = c(0.5, 1.3, 1.9),
                           plot_seg = TRUE)
```

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
Obj_filtered$ref=Obj_filtered$select_normal$region_normal[8] # choose one normal region
```
<br/>

#### Step5. Genotype each cell in each region

* Estimate cell-specific (rho_hat, theta_hat) values for each region.
```
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, ref_counts = ref_counts) # for cell line without normal cells in the tumor sample.
```

* Genotype all cells and generate a genotype plot for each region.
```
Obj_filtered=Genotype(Obj_filtered = Obj_filtered)
```
The output genotying results for the five regions are shown below.

<br/><br/>

#### Step6. Construct lineage structure using genotypes for each cell across all regions

* Generate lineage tree based on cell-specific genotypes across the regions.
```
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 10)
```
The output clustering result for the five regions is shown below.

<br/>

## Citation
Wu, C.-Y. et al. Alleloscope: Integrative analysis of single cell haplotype-divergent copy number alterations and chromatin accessibility changes reveals novel clonal architecture of cancers. bioRxiv (2020): [https://doi.org/10.1101/2020.10.23.349407](https://doi.org/10.1101/2020.10.23.349407)

