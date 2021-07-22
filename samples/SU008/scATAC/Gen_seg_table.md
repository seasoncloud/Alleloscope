Alleloscope (WES/WGS)
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description
Alleloscope provides a function to generate the segmentation table (seg_table) for the downstream analysis from matched WGS or WES data. 

For more information about the method, please check out the [github](https://github.com/seasoncloud/Alleloscope) and the [paper](https://doi.org/10.1038/s41587-021-00911-w).
<br/>

## Prepare input files
Two BED files (for tumor and normal WGS/WES) with each row a genomic bin. The first column is chromosomes; the second column is start sites; thee third column is end sites; the fourth column is summed read counts. 
This can be generated using bedtools with the following command.
* For GRCh38 reference genome and 100kb bins, the 'hg38.100Kb.windows.sorted.bed' file [here](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/hg38.100Kb.windows.sorted.bed) can be used. 
```
# for tumor sample
bedtools intersect -a ~/hg38.100Kb.windows.sorted.bed \
-b ~/sample.bam \
-c -sorted \
> ~/hg38.100Kb.windows.counts.tumor.bedg

# for germline sample
bedtools intersect -a ~/hg38.100Kb.windows.sorted.bed \
-b ~/sample.bam \
-c -sorted \
> ~/hg38.100Kb.windows.counts.germline.bedg
```

#### Step0. Load the input files

* In R, set up the environment and read files
```
library(Alleloscope) # load the library
setwd("~/Alleloscope/") # set path to the github folder

dir_path <- "./WGS/"; dir.create(dir_path) # set up output directory

size=read.table("data-raw/sizes.cellranger-atac-hg19-1.2.0.txt", stringsAsFactors = F) # read size file
size=size[1:22,]

WGSt=read.table("~/hg38.100Kb.windows.counts.tumor.bedg", stringsAsFactors = F, sep='\t')
WGSn=read.table("~/hg38.100Kb.windows.counts.germline.bedg", stringsAsFactors = F, sep='\t')
WGSt=WGSt[which(WGSt$V1 %in% paste0('chr',1:22)), ]
WGSn=WGSn[which(WGSn$V1 %in% paste0('chr',1:22)), ]
WGSt[,2]=as.numeric(WGSt[,2])+1; WGSt[,3]=as.numeric(WGSt[,3])+1
WGSn[,2]=as.numeric(WGSn[,2])+1; WGSn[,3]=as.numeric(WGSn[,3])+1
rownames(WGSn)=rownames(WGSt)=paste0(WGSt$V1,"-", WGSt$V2,'-', WGSt$V3)
WGSt=WGSt[,4, drop=F]; WGSn=WGSn[,4, drop=F]
```

#### Step1. Creating a Alleloscope object for the analysis

* First, create an Alleloscope obj
```
Obj_filtered=Createobj(samplename=name, genome_assembly="GRCh38", dir_path=dir_path, size=size, assay='WGS', cov=TRUE)
```
<br/>

#### Step2. Unbiased segmentation based on matched WES/WGS data

* Load the segmentation results of the matched WES data from Yost et al., 2019.
```
Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                          raw_counts=WGSt, # from matched DNA sequencing (bulk/single)
                          ref_counts=WGSn, # from matched DNA sequencing (bulk/single)
                          plot_seg = TRUE)
```
<br/>
* The segmentation plot and seg_table.rds is stored in the dir_path.


## Citation
Wu, C.-Y. et al. Integrative single-cell analysis of allele-specific copy number alterations and chromatin accessibility in cancer. Nature Biotechnology (2021): [https://doi.org/10.1038/s41587-021-00911-w](https://doi.org/10.1038/s41587-021-00911-w)





[Back to the main page](https://github.com/seasoncloud/Alleloscope)

