#' Select the segments in the "seg_table" with more than nSNP
#'
#' @param Obj_filtered An Alleloscope object with SNP info and raw segmentation table "seg_table".
#' @param nSNP An integer of minimum number of SNPs for region selecgtion.
#' @param fvar Logical(TRUE/FALSE) Whether or not to filter out segments with variance >0.9 quantile. 
#'
#' @return A Alleloscope object with "seg_table_filtered" added.
#'
#' @export
Segments_filter=function(Obj_filtered=NULL, nSNP=0 ,len=500000, fvar=FALSE){
  assay=Obj_filtered$assay
  seg_table=Obj_filtered$seg_table
  var_list=Obj_filtered$var_all

  var_str=paste0(as.character(var_list[,1]),":", as.character(var_list[,2]),"_", as.character(var_list[,4]),"_", as.character(var_list[,5]))
  var_chr=as.numeric(sapply(strsplit(var_str,':'),'[',1))
  var_pos=as.numeric(sapply(strsplit(sapply(strsplit(var_str,':'),'[',2), "_"),'[',1))

  seg_table=data.frame(seg_table, stringsAsFactors = F)
  query=GenomicRanges::GRanges(seg_table$chr, IRanges::IRanges(as.numeric(seg_table$start)+1,as.numeric(seg_table$end))) ## cytoband 0-based start and 1-based end
  subject=GenomicRanges::GRanges(var_chr, IRanges::IRanges(var_pos, var_pos))
  ov=findOverlaps(query, subject)
  ov=as.matrix(ov)
  seg_table=cbind(seg_table, table(ov[,1])[match(1:(dim(seg_table)[1]), names(table(ov[,1])))])
  row.names(seg_table)=NULL
  chrr= paste0(seg_table$chr,":", seg_table$start)
  seg_table=cbind(seg_table,chrr)
  seg_table=seg_table[!is.na(seg_table$Freq),]
  seg_table$chrr=as.character(seg_table$chrr)

  if(nrow(seg_table)>22 & fvar==TRUE){
  ind=as.numeric(seg_table$Var1[which(as.numeric(seg_table$Freq)>nSNP & as.numeric(seg_table$length)>len & as.numeric(seg_table$var)<quantile(as.numeric(seg_table$var), 0.99,na.rm = T))])
  }else{
    ind=as.numeric(seg_table$Var1[which(as.numeric(seg_table$Freq)>nSNP &  as.numeric(seg_table$length)>len)])
  }
  seg_table=seg_table[ind,]

  Obj_filtered[['seg_table_filtered']]=seg_table
  Obj_filtered[['nSNP']]=nSNP

  message("Segments have been filtered!")
  cat("\"seg_table_filtered\" was added to the Obj_filtered object.\n")
  return(Obj_filtered)

}
