#' Generate Alleloscope object for analysis
#'
#' @param alt_all A SNP by cell read count matrix/ spare matrix for the alternative alleles.
#' @param ref_all A SNP by cell read count matrix/ spare matrix for the reference alleles.
#' @param vcf_all A matrix/ data.frame of the vcf format for SNP information. (The length and order are the same as nrow(alt_all) and nrow(alt_all))
#' @param samplename Sample name for the data.
#' @param genome_assembly The genome assembly used for sequencing alignment. (ex: "GRCh38" or "GRCh37")
#' @param dir_path Path of the output directory.
#' @param barcodes A matrix/ data.frame with barcodes for each cell in the first column.
#' @param size A matrix with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes. 
#' @param assay A character indicating the type of sequencing data. (ex: "scDNAseq" or "scATACseq")
#' @param cov Logical (TRUE/FALSE). Whether or not to use only coverage. If "cov" is TRUE, alt_all, ref_all, andvcf_all are not required.
#' 
#' @import Matrix
#' @return A Alleloscope object including the necessary information.
#'
#' @export
Createobj=function(alt_all=NULL, ref_all=NULL, var_all=NULL, samplename='sample',genome_assembly="GRCh38", dir_path='./', barcodes=NULL, size=NULL, assay='scDNAseq', cov=FALSE){
  
  # check parameters
 if(!(nrow(size)>0 & ncol(size)==2)){
    stop("Please provide a matrix/ data.frame with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.")
  }
  
  dir.create(paste0(dir_path,"/plots"))
  
  
  if(cov==FALSE){
    if(!(nrow(barcodes)>0 & ncol(barcodes)==1)){
      stop("Please provide a matrix/ data.frame with barcodes for each cell in the first column.")
    }
    
    if(!grepl('chr',var_all[1,1])){
      var_all[,1]=paste0('chr', var_all[,1])
    }
    
    ##
    colnames(alt_all)=barcodes[,1]
    colnames(ref_all)=barcodes[,1]
    
    total_all=alt_all+ref_all
    #dim(total_all)
  
  ## read the size file
  
  if(grepl('chr', as.character(size[2,1]))){
    size[,1]=sapply(strsplit(size[,1],'hr'),'[',2)##if chr
  }
  size=size[which(size[,1] %in% as.character(1:22)),]
  size=size[order(as.numeric(as.character(size$V1))),]
  size_name=size[,1]
  size=as.numeric(as.character(size[,2]))
  names(size)=size_name
  
  size=size[which(paste0('chr',size_name) %in% unique(var_all$V1))]
  
  ## read the meta data
  #cell_info=read.table(path_cell_summary, sep=',', header = T, stringsAsFactors = F)
  
  
  
  output=list("alt_all"=alt_all, "ref_all"=ref_all, "total_all"=total_all,"var_all"=var_all, "barcodes"=barcodes[,1], "size"=size,"samplename"=samplename,
              "dir_path"=dir_path,"genome_assembly"=genome_assembly,
              #"cell_info" = cell_info,
              "cell_filter"=NULL, "SNP_filter"=NULL, "min_vaf"=NULL, "max_vaf"=NULL,
              "seg_table"=NULL, "seg_table_filtered"=NULL, "nSNP"=NULL, "rds_list"=NULL, "select_normal"=NULL,
              "ref"=NULL, "genotype_table"=NULL, "assay"=assay)
  
}else{
  if(grepl('chr', as.character(size[2,1]))){
    size[,1]=sapply(strsplit(size[,1],'hr'),'[',2)##if chr
  }
  size=size[which(size[,1] %in% as.character(1:22)),]
  size=size[order(as.numeric(as.character(size$V1))),]
  size_name=size[,1]
  size=as.numeric(as.character(size[,2]))
  names(size)=size_name
  
  output=list("size"=size,"samplename"=samplename,
              "dir_path"=dir_path,"genome_assembly"=genome_assembly,
              #"cell_info" = cell_info,
              "cell_filter"=NULL, "SNP_filter"=NULL, "min_vaf"=NULL, "max_vaf"=NULL,
              "seg_table"=NULL, "seg_table_filtered"=NULL, "nSNP"=NULL, "rds_list"=NULL, "select_normal"=NULL,
              "ref"=NULL, "genotype_table"=NULL, "assay"=assay)
  
}
  
  
  message("Object successfully created!")
  return(output)
  
}
