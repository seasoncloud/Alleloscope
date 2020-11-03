#' Identify candidate normal cells and normal regions for cell coverage normalization
#'
#' @param Obj_filtered An Alleloscope object with major haplotype proportion (theta_hat) for each cell of each region in the "rds_list".
#' @param raw_counts A large binned coverage matrix (bin by cell) with values being read counts for all chromosomal regions of tumor sample.
#' @param cell_nclust Integer. Number of clusters used in identifying normal cells in the sample.
#' @param plot_theta Logical (TRUE/FALSE). Whether or not to plot the hierarchical clustering result using the theta_hat values across regions.
#' @param cell_type A matrix with two columns: COL1- cell barcodes; COL2- cell types ("tumor" and others)
#' @param cutree_rows Integer. Number of clusters the rows are divided into for visualization (inherited from the pheatmap function).
#'
#' @return A Alleloscope object with a "select_normal" list added.
#' A "select_normal" list includes
#' "barcode_normal": Barcodes of the identified normal cells in the tumor sample.
#' "region_normal": A vector of ordered potential normal regions for selection. (1st is the most possible.)
#' "region_normal_rank": A table with the potential "normal regions" for the k clusters from hierarchical clustering.
#' "k_normal": An integer indicates the kth clsuter that is idenfied as "normal cells"
#'
#' @export
Select_normal=function(Obj_filtered=NULL, raw_counts=NULL, cell_nclust=5 , plot_theta=FALSE, cell_type=NULL, cutree_rows=3 ){

EMresult=Obj_filtered$rds_list
filtered_seg_table=Obj_filtered$seg_table_filtered
#cell_info=Obj_filtered$cell_info

plot_path=paste0(Obj_filtered$dir_path,'/plots/')

samplename=Obj_filtered$samplename
cell_barcodes=Obj_filtered$barcodes
assay=Obj_filtered$assay
ncell=length(Obj_filtered$barcodes)
distance_seg=as.numeric(filtered_seg_table$end)-as.numeric(filtered_seg_table$start) ## for cytoarm
size=Obj_filtered$size
cell_total=Matrix::colSums(Obj_filtered$total_all)

## raw/ref count matrix info
raw_chr=sapply(strsplit(rownames(raw_counts),'-'),'[',1)
raw_start=as.numeric(sapply(strsplit(rownames(raw_counts),'-'),'[',2))
raw_end=as.numeric(sapply(strsplit(rownames(raw_counts),'-'),'[',3))

theta_N=list()

for(chrr in as.character(filtered_seg_table$chrr)){ # for cytoarm

  chrrn=unlist(strsplit(chrr,':'))[1]
  result=EMresult[[paste0('chr', chrr)]]
  theta_hat=result$theta_hat
  names(theta_hat)=result$barcodes
  barcodes=result$barcodes#

  # subset the coverge info
  raw_counts_chr=raw_counts[which(raw_chr %in% paste0('chr', as.character(chrrn))),]
  raw_chr_sub=raw_chr[which(raw_chr %in% paste0('chr', as.character(chrrn)))]
  raw_start_sub=raw_start[which(raw_chr %in% paste0('chr', as.character(chrrn)))]
  raw_end_sub=raw_end[which(raw_chr %in% paste0('chr', as.character(chrrn)))]
  
  query=GenomicRanges::GRanges(paste0('chr',chrrn),IRanges::IRanges(as.numeric(filtered_seg_table[which(filtered_seg_table$chrr == chrr),2]), as.numeric(filtered_seg_table[which(filtered_seg_table$chrr == chrr),3])))
  subject=GenomicRanges::GRanges(raw_chr_sub, IRanges::IRanges(as.numeric(raw_start_sub),as.numeric(raw_end_sub))) ## cytoband 1-based start and 1-based end
  ov=findOverlaps(query, subject)
  ov=as.matrix(ov)
  
  
  bin_start=min(ov[,2])
  bin_end=max(ov[,2])
  raw_counts_chr=raw_counts_chr[bin_start:bin_end,] ## for p and q #for cytoarm

  raw_counts_chr=raw_counts_chr[, 1:length(Obj_filtered$barcodes)][,match(barcodes, cell_barcodes)]
  Nr=colSums(raw_counts_chr)

 # N0=cell_info$num_mapped_dedup_reads
  N0=cell_total

  barcodes_non_noisy=cell_barcodes#cell_info$barcode[which(cell_info$is_noisy==0)]
  barcodes_non_noisy=intersect(barcodes_non_noisy, barcodes)

  names(Nr)=barcodes
  Nr=Nr[match(barcodes_non_noisy, barcodes)]
  N0=N0[match(barcodes_non_noisy, barcodes)]
  Ni=Nr/N0*sum(size)/distance_seg[which(filtered_seg_table$chrr==chrr)]
  names(Ni)=barcodes_non_noisy
  Ni[Ni>quantile(Ni, 0.99)]=quantile(Ni, 0.99)


  theta_hat=theta_hat[match(barcodes_non_noisy, names(theta_hat))]


  df=data.frame(Nir_Ni0=Ni, theta_hat=theta_hat)

  theta_N[[paste0("rho_",as.character(chrr))]]=Ni
  theta_N[[paste0("theta_",as.character(chrr))]]=theta_hat

  cat(paste0(chrr," "))
}

cat("\n")
#####3 remove the region that have too few cells!

cell_list<-lapply(theta_N, function(x) {
  names(x)
})

cell_intersect <- Reduce(intersect, cell_list)


theta_hat_cbn <- sapply(theta_N,function(x){
  x[match(cell_intersect, names(x))]
})

rownames(theta_hat_cbn) <- cell_intersect
saveRDS(theta_hat_cbn, paste0(Obj_filtered$dir_path,"/rds/theta_N_seg.rds"))


theta_hat_cbn2=theta_hat_cbn[,which(stringr::str_sub(colnames(theta_hat_cbn), end=1)=='t'), drop=F]
rho_hat_cbn2=theta_hat_cbn[,which(stringr::str_sub(colnames(theta_hat_cbn), end=1)=='r'), drop=F]

tmp=pheatmap::pheatmap(theta_hat_cbn2, cluster_cols = F, cluster_rows = T, show_rownames = F, clustering_method = "ward.D2", silent = TRUE)#, gaps_col=gaps_col), annotation_row = cell_label, annotation_col =region_label)


if(plot_theta==TRUE){
pdf(paste0(plot_path,"/hierarchcial_clustering_theta.pdf"), width = 12,height = 6)
if(ncol(theta_hat_cbn2)<2){
  gaps_col=NULL
}else{
  gaps_col=2*(1:(dim(theta_hat_cbn2)[2]/2))
}
  region_name=paste0('chr',sapply(strsplit(colnames(theta_hat_cbn2),"_"),'[',2))
if(is.null(cell_type)){
tmp=pheatmap::pheatmap(theta_hat_cbn2, cluster_cols = F, cluster_rows = T, show_rownames = F, clustering_method = "ward.D2", gaps_col=gaps_col,  labels_col=region_name)#, annotation_row = cell_label, annotation_col =region_label)
}else{
  barcodes_tumor=cell_type[which(cell_type[,2]=='tumor'),1]
  barcodes_normal=cell_type[which(cell_type[,2]!='tumor'),1]
  cell_label=rep('cell', dim(theta_hat_cbn2)[1])
  cell_label[which(rownames(theta_hat_cbn2) %in% barcodes_tumor)]='tumor'
  cell_label[which(rownames(theta_hat_cbn2) %in% barcodes_normal)]='normal'
  cell_label=data.frame('cell type'=cell_label, row.names = rownames(theta_hat_cbn2))
  my_colour = list(
    cell.type = c(tumor = "#d53e4f", normal="#1f78b4")
  )
  
tmp=pheatmap::pheatmap(theta_hat_cbn2, cluster_cols = F, cluster_rows = T, show_rownames = F, clustering_method = "ward.D2", gaps_col=1:(dim(theta_hat_cbn2)[2]), annotation_row = cell_label, annotation_colors = my_colour, cutree_rows = cutree_rows,  labels_col=region_name)
}
dev.off()

cat(paste0("Plot for theta_hat clustering across all regions is saved in the path:", plot_path,"\n"))
}

## select normal cells
k=cell_nclust
clust=cutree(tmp$tree_row, k=k)
theta_ss_k=c()
for(ii in 1:k){
  theta_sub=theta_hat_cbn2[which(clust==ii),, drop=F]
  theta_mean=apply(theta_sub, 2, mean)
  theta_ss=sum((theta_mean-0.5)^2)
  theta_ss_k=c(theta_ss_k, theta_ss)
}

(k_normal=which(theta_ss_k==min(theta_ss_k)))

barcode_normal=names(clust)[which(clust==k_normal)]




## select normal regions
#k=5)
region_normal_rank5=matrix(nrow = 10, ncol=k)
for(ii in 1:k){
  theta_sub=theta_hat_cbn2[which(clust==ii),, drop=F]
  rho_sub=rho_hat_cbn2[which(clust==ii),, drop=F]
  theta_ss_region=apply(theta_sub, 2, function(x) sum((x-0.5)^2))
  names(theta_ss_region)=colnames(theta_hat_cbn2)
  #rho_cv_region=apply(rho_hat_cbn2, 2, function(x) sd(x)/mean(x))
  rho_med_region=apply(rho_sub,2, mean)
  normal_rank=sort(rank(theta_ss_region)+rank(rho_med_region))
  normal_regions=sapply(strsplit(names(normal_rank),"_"),'[',2)
  region_normal_rank5[,ii]=normal_regions[1:10]
}


colnames(region_normal_rank5)=paste0('k',1:k)
rownames(region_normal_rank5)=paste0('rank',1:10)

tmp=c()
for(rr in 1:nrow(region_normal_rank5)){
  tmp=c(tmp, names(table(region_normal_rank5[rr,-k_normal])))
}

(region_normal=names(sort(table(tmp), decreasing = T))[1:10])

select_normal=list("barcode_normal"=barcode_normal, "region_normal"=region_normal, "region_normal_rank"=region_normal_rank5, "k_normal"=k_normal )
Obj_filtered$select_normal=select_normal
message("Candidate normal cell and normal region info is in \"Obj_filtered$select_normal\".")
return(Obj_filtered)

}
