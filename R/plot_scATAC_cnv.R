#' Generate Alleloscope object for analysis
#'
#' @param raw_mat A binned coverage matrix (m1 bin by n1 cell) with values being read counts for scATAC-seq of tumor sample (with some normal cells). The matrix can be generated using https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/Gen_bin_cell_atac.R
#' @param cell_type A matrix with two columns: COL1- cell barcodes; COL2- cell types (Tumor cells should be labeled with "tumor1, tumor2 and etc.").
#' @param normal_lab Character(s) indicating the cell types considered as normal cells. If not specify, "normal" cell type should exist in the cell_type dataframe. 
#' @param size A matrix with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes (eq.1-22). 
#' @param window_w window size for signal pooling in individual cells.
#' @param window_step step size for signal smoothing in individual cells.
#' @param plot_path Path to plot the heatmap.
#' @param nclust Integer. Number of clusters the rows are divided into for visualization and clustering.
#' @param var.filter Logical (TRUE/FALSE) Whether or not to filter our highly variable features.
#' 
#' @import matrixStats
#' @return A vector indicating the ordered cluster number (from hierarchical clustering) of each cell and a heatmap saved.
#'
#' @export
plot_scATAC_cnv=function(raw_mat=NULL,cell_type=NULL,normal_lab="normal", size=NULL, window_w=10000000, window_step=2000000, plot_path=NULL, nclust=3, var.filter=FALSE){
  
  if(is.null(plot_path)){
    plot_path=paste0("./plots/CNV_cov_w",window_w,"_s",window_step,"_sub.pdf")
    dir.create("./plots/")
  }
  ## normlize by cell size
  if(var.filter==TRUE){
  vars=apply(raw_mat,1, var)
  cellsize=colSums(raw_mat[which(vars<quantile(vars,0.99)),])
  }else{
    cellsize=colSums(raw_mat)
  }
  cellsize=matrix(rep(cellsize, nrow(raw_mat)), byrow =T, ncol=ncol(raw_mat))
  raw_mat=raw_mat/cellsize
  if(grepl("chr",size[1,1])){
    size=size
  }else{
    size[,1]=paste0("chr",size[,1])
  }
  size=size[ (size[,1] %in% paste0('chr',1:22)),]
  cnv_bin0=GRanges(size[,1], IRanges(1,as.numeric(size[,2])))
  sw=slidingWindows(x = cnv_bin0, width = window_w, step = window_step)
  #width(sw)
  cnv_bin=sw@unlistData
  
  subject=GRanges(sapply(strsplit(rownames(raw_mat),":|_|-"),'[',1),
                  IRanges(as.numeric(sapply(strsplit(rownames(raw_mat),":|_|-"),'[',2))+1,
                          as.numeric(sapply(strsplit(rownames(raw_mat),":|_|-"),'[',3))))
  ov=findOverlaps(cnv_bin, subject )
  ov=as.matrix(ov)
  
  smooth_mat=matrix(ncol=ncol(raw_mat), nrow=length(cnv_bin))
  for(ii in 1:length(cnv_bin)){
    ri=colSums(raw_mat[ov[which(ov[,1]==ii),2],, drop=F])
    smooth_mat[ii,]=ri
  }
  colnames(smooth_mat)=colnames(raw_mat)
  rownames(smooth_mat)=paste0(as.character(seqnames(cnv_bin)),'-',start(cnv_bin),'-', end(cnv_bin))
  rownames(cell_type)=cell_type[,1]
  mat_celltype=cell_type[match(colnames(smooth_mat), rownames(cell_type)),2]
  #tumor_mat=smooth_mat[,which(!grepl('tumor',mat_celltype))]###
  nontumor_mat=smooth_mat[,which((mat_celltype %in% normal_lab))]###
  #tumor_mat=smooth_mat[,which(mat_celltype=='tumor')]
  #nontumor_mat=smooth_mat[,which(mat_celltype!='tumor')]
  
  
  bin_ind=which(rowMedians(nontumor_mat)!=0)
  norm_matrix=matrix(rep(rowMedians(nontumor_mat),ncol(smooth_mat)), ncol=ncol(smooth_mat), byrow=F)
  smooth_mat_norm=smooth_mat[bin_ind,]/(norm_matrix[bin_ind,])
  smooth_mat_norm=apply(smooth_mat_norm, c(1,2), function(x) min(x,5))
  plot_matrix=t(smooth_mat_norm)
  chrgap=(table(sapply(strsplit(rownames(smooth_mat[bin_ind,]),'-'),'[',1))[paste0('chr',1:nrow(size))])
  chrgap[is.na(chrgap)]=0
  col_lab=rep(" ", ncol(plot_matrix))
  #col_lab[c(0, cumsum(chrgap)[1:(length(chrgap)-1)])+chrgap/2]=paste0("chr",as.character(1:nrow(size)))
  col_lab[c(0, cumsum(chrgap)[which(chrgap!=0)])[1:length(which(chrgap!=0))]+chrgap[which(chrgap!=0)]/2]=names(chrgap)[!is.na(names(chrgap))]
  celltype=cell_type
  #rownames(celltype)=celltype[,1]
  #celltype=celltype[,-1, drop=F]
  celltype=as.data.frame(cell_type)
  rownames(celltype)=celltype[,1]
  
  ## break for coloring
  plot_matrix=plot_matrix*2
  plot_matrix=apply(plot_matrix, c(1,2), function(x) if(x<=2.5 & x>=1.5){x=2}else{x=x})
  
  breaklength = 50
  setcolor = colorRampPalette(c("blue", "white", "red"))(breaklength)
  setbreaks = c(seq(min(plot_matrix), 1.7, length.out=ceiling(breaklength/2) + 1), 
                c(2.3,seq((max(plot_matrix)-2.3)/breaklength+2.3, max(plot_matrix), 
                          length.out=floor(breaklength/2)))[1:(breaklength/2)])
  
  
  #pdf(plot_path, width=16, height=9)
  png(plot_path,width=800, height=450)
  tmp=pheatmap::pheatmap(plot_matrix,
                         cluster_cols = F, cluster_rows = TRUE,
                         show_rownames = F,
                         show_colnames = T,
                         color=setcolor,
                         breaks = setbreaks,
                         labels_col=col_lab,
                         clustering_distance_rows = "correlation",
                         clustering_method = "ward.D2",
                         gaps_col=cumsum(chrgap),
                         cutree_rows = nclust,
                         annotation_row=celltype[, ncol(celltype),drop=F])
  
  dev.off()
  
  message("Plot successfully generated!")
  cat("Path: ",plot_path,"\n")
  
  od=tmp$tree_row$order
  clust=cutree(tmp$tree_row, k=nclust)
  clust_order=clust[od]
  
  
  cov_obj=list(clust_order=clust_order, plot_matrix=plot_matrix)
  return(cov_obj)
  
  
}
