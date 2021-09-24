#' Run all steps for scDNA-seq data
#'
#' @param alt_all A SNP by cell read count matrix/ spare matrix for the alternative alleles.
#' @param ref_all A SNP by cell read count matrix/ spare matrix for the reference alleles.
#' @param vcf_all A matrix/ data.frame of the vcf format for SNP information. (The length and order are the same as nrow(alt_all) and nrow(alt_all))
#' @param samplename Sample name for the data.
#' @param genome_assembly The genome assembly used for sequencing alignment. (ex: "GRCh38" or "GRCh37")
#' @param dir_path Path of the output directory.
#' @param barcodes A matrix/ data.frame with barcodes for each cell in the first column.
#' @param size A numeric vector for the size (bp) of different chromosomes (with the names indicating which chromosome from 1 to 22)
#' @param assay A character indicating the type of sequencing data. (ex: "scDNAseq" or "scATACseq")
#' @param cell_filter An integer of minimum cell number for SNP selection.
#' @param SNP_filter An integer of minimum SNP number for cell selection.
#' @param min_vaf A numerical value in the range (0,1) of minimum SNP variant allele frequency in the pseudo bulk for SNP selection.
#' @param max_vaf A numerical value in the range (0,1) of mzsimum SNP variant allele frequency in the pseudo bulk for SNP selection.
#' @param raw_counts A large binned coverage matrix (m1 bin by n1 cell) for all chromosomal regions of tumor sample.
#' @param ref_counts A large binned coverage matrix (m2 bin by n2 cell) for all chromosomal regions of normal sample.
#' @param type Specify whethere the sample is a "tumor" or "cellline". If "type" is a "cellline", param "ref_counts" needs to be specified for normal sample.
#' 
#' @import Matrix
#' @return A Alleloscope object including the necessary information.
#'
#' @export
Rundf_dna=function(alt_all=NULL, ref_all=NULL, var_all=NULL, samplename='sample',genome_assembly="GRCh38", 
                   dir_path='./', barcodes=NULL, size=NULL, assay='scDNAseq',
                   raw_counts=NULL, ref_counts=NULL, type='tumor', 
                   cell_filter=5, SNP_filter=10 ,min_vaf=0, max_vaf=1){
  
  Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename=samplename, genome_assembly=genome_assembly, dir_path=dir_path, barcodes=barcodes, size=size, assay=assay)
  
  Obj_filtered=Matrix_filter(Obj=Obj, cell_filter = cell_filter, SNP_filter = SNP_filter, min_vaf = min_vaf, max_vaf = max_vaf)
  
  rm(Obj)
  
  Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                            raw_counts=raw_counts,
                            ref_counts=ref_counts,
                            plot_seg = TRUE)
  
  Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP = 200)
  
  Obj_filtered=Est_regions(Obj_filtered = Obj_filtered)
  
  Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts, plot_theta = T)
  
  Obj_filtered$ref=Obj_filtered$select_normal$region_normal[1] 
  
  Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type=type, raw_counts=raw_counts, ref_counts = ref_counts,cov_adj=1) 
  
  Obj_filtered=Genotype(Obj_filtered = Obj_filtered)
  
  linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 5)
  
  saveRDS(Obj_filtered,paste0(dir_path, "rds/Obj_filtered.rds"))
  
  return(Obj_filtered)
  message("Run has been completed!")
  
}
