#' Perform EM iterations on the filtered cells with barcodes, and plot the results for each region.
#'
#' @param Obj_filtered An Alleloscope object with allele and segment information for estimating cell major haplotype proportion (theta_hat) for each region.
#' @param max_nSNP Integer. Maximum SNP number used for estimating theta_hat for a region.
#' @param plot_stat Logical (TRUE/ FALSE). Whether or not to plot the statistics and EM results for each region.
#' @param min_ncell Integer. Filter out the cells with reads < min_ncells.
#' @param rds_path The path for saving the rds files for the estimated results for each region.
#' @param cont Logical (TRUE/FALSE). Whether or not to skip the regions with the rds files in the specified rds_path.
#' @param max_iter Integer. Maximum numbers of EM iterations. 
#' @param phases List. The estimated phase indicators (I_j) of each SNP across all regions.
#'
#' @return A "rds_list" of the estimated SNP phases (I_hat), estimated cell major haplotype proportion (theta_hat) for all regions.
#'
#' @import rtracklayer
#' @export
Est_regions=function(Obj_filtered=NULL,max_nSNP=30000, plot_stat=TRUE,min_ncell=20, rds_path=NULL, cont=FALSE, max_iter=50, phases=NULL){
  message("Estimation starts.")
  assay=Obj_filtered$assay
  filtered_seg_table=Obj_filtered$seg_table_filtered
  #cell_info=Obj_filtered$cell_info

  plot_path=paste0(Obj_filtered$dir_path,'/plots/')
  samplename=Obj_filtered$samplename
  if(is.null(rds_path)){
  rds_path=paste0(Obj_filtered$dir_path,"/rds")
  }
  dir.create(rds_path)
  
  dir.create(paste0(plot_path, "EMresults"))
  dir.create(paste0(rds_path, "/EMresults"))


  if(is.null(filtered_seg_table)){
    message("Estimation for each chromosome")
    filtered_seg_table=data.frame("chr"=1:length(size), 'start'=rep(0,length(size)), 'end'=Obj_filtered$size, 'chrr'=as.character(1:length(size)), stringsAsFactors = F)
  }

  ## look at segments
  var_list = Obj_filtered$var_all

  var_str=paste0(as.character(var_list[,1]),":", as.character(var_list[,2]),"_", as.character(var_list[,4]),"_", as.character(var_list[,5]))
  var_chr=as.numeric(sapply(strsplit(var_str,':'),'[',1))
  var_pos=as.numeric(sapply(strsplit(sapply(strsplit(var_str,':'),'[',2), "_"),'[',1))



  query=GRanges(filtered_seg_table$chr, IRanges(as.numeric(filtered_seg_table$start)+1,as.numeric(filtered_seg_table$end))) ## cytoband 0-based start and 1-based end
  subject=GRanges(var_chr, IRanges(var_pos, var_pos))
  ov=findOverlaps(query, subject)
  ov=as.matrix(ov)

  selseg=filtered_seg_table$chrr


 # if(non_noisy==TRUE){
#    barcodes_non_noisy=cell_info$barcode[which(cell_info$is_noisy==0)]}

  rds_list=list()


  for(chrr in as.character(selseg)){
    
    if(cont==FALSE){
      chk="chk"
    }else{
      chk=paste0('chr',chrr,'.rds')
    }

    if(! chk %in% list.files(rds_path)){

    chr_ind=ov[which(ov[,1] %in% which(filtered_seg_table$chrr %in% chrr)),2]

    alt_all_sub=Obj_filtered$alt_all[chr_ind,, drop=F]
    total_all_sub=Obj_filtered$total_all[chr_ind,, drop=F]
    var_list_sub=var_list[chr_ind,]

   # if(non_noisy==TRUE){
  #    alt_all_sub=alt_all_sub[,which(colnames(alt_all_sub) %in% barcodes_non_noisy), drop=F]
   #   total_all_sub=total_all_sub[,which(colnames(total_all_sub) %in% barcodes_non_noisy), drop=F]
  #  }


    cc=Matrix::colSums(total_all_sub)
    cc_ind_sub=which(cc>min_ncell)
    
    alt_all_sub=alt_all_sub[,cc_ind_sub,drop=F]
    total_all_sub=total_all_sub[,cc_ind_sub,drop=F]
    #var_list_sub=var_list_sub[rr_ind_sub,1:5]
    
    if(dim(total_all_sub)[2]==0){
      cat(paste0("No cells after filtering for ", chrr, " "))
      next
    }
    
    
    rr_ind_sub=1:nrow(total_all_sub)
    rm_ind=which(Matrix::rowSums(total_all_sub)==0)
    if(length(rm_ind)!=0){
    rr_ind_sub=rr_ind_sub[-rm_ind]
    }
    
    if(length(rr_ind_sub)>max_nSNP){
      rr_ind_sub=sort(sample(rr_ind_sub,max_nSNP))}
    
    alt_all_sub=alt_all_sub[rr_ind_sub, ,drop=F]
    total_all_sub=total_all_sub[rr_ind_sub, ,drop=F]
    var_list_sub=var_list_sub[rr_ind_sub,1:5]
    
    
    
    af_all_sub=Matrix::rowSums(alt_all_sub)/Matrix::rowSums(total_all_sub)
    
    af_all_sub[is.na(af_all_sub)]=0
    
    if(dim(alt_all_sub)[1]<3 | length(unique(af_all_sub))<3){
      cat(paste0("No SNPs after filtering for ", chrr, " "))
      next
    }

 ## plot stats for each region
    
    if(plot_stat==TRUE){
      pdf(paste0(plot_path,"/EMresults/statistics_",assay,"_chr", as.character(chrr), '.pdf' ))
      par(mfrow=c(3,1))
      hist(Matrix::colSums(total_all_sub), main=paste0(samplename," ",assay," chr",as.character(chrr)," coverage (",as.character(dim(total_all_sub)[2])," cells across ",  as.character(dim(total_all_sub)[1]), " SNPs)"), xlab="coverage of individul cells", ylab='frequency', breaks = 100)
      hist(Matrix::rowSums(total_all_sub), main=paste0(samplename," ",assay," chr",as.character(chrr)," coverage (",as.character(dim(total_all_sub)[1])," SNPs across ",  as.character(dim(total_all_sub)[2]), " cells)"), xlab="coverage of individul SNPs", xlim=c(0,100), ylab='frequency', breaks = 1000)
      hist(af_all_sub, 100, main="Histogram of VAF values")
      dev.off()
    }

  if(is.null(phases)){
    #message("EM iterations for each region.")
    result=EM(ref_table = as.matrix(total_all_sub-alt_all_sub), alt_table = as.matrix(alt_all_sub) ,seed=1000, max_iter=max_iter)

    result$barcodes=colnames(total_all_sub)
    result$SNPs=paste0('chr', var_list_sub$V1,':', var_list_sub$V2,'_', var_list_sub$V4,'_', var_list_sub$V5)
  }else{
    #message("Estimate theta_i using phasing info.")
    result=list()
    #dna_snp=unlist(Obj_scDNA$rds_list[[paste0('chr',chrr)]][5])
    dna_ind=unlist(phases[[paste0('chr',chrr)]][['I_hat']])
    #dna_ind=as.numeric(phases)
    atac_snp=paste0('chr',var_list_sub$V1,':', var_list_sub$V2,'_', var_list_sub$V4,'_', var_list_sub$V5)
    #inter=intersect(dna_snp, atac_snp)
    
  #  alt_all_sub=alt_all_sub[match(inter, atac_snp),]
  #  total_all_sub=total_all_sub[match(inter, atac_snp),]
  #  var_list_sub=var_list_sub[match(inter, atac_snp),]
    
  #  ind=as.numeric(dna_ind[match(inter, dna_snp)])
    ind=dna_ind
    
    alt_table=t(as.matrix(alt_all_sub))
    tot_table=t(as.matrix(total_all_sub))
    ref_table=tot_table- alt_table  
    mm=dim(alt_table)[2]  # mm snv
    nn=dim(alt_table)[1]  # nn cell
    
    ind_table=matrix(rep(ind, nn), nrow = nn, byrow = T)
    
    # maximization step
    w1=rowSums((ref_table*ind_table)+(alt_table*(1-ind_table)))
    w2=rowSums((ref_table*(1-ind_table))+(alt_table*ind_table))
    
    theta=w1/(w1+w2)
    theta[is.na(theta)]=0.5
    result=list("theta_hat"=theta, "I_hat"=ind, "barcodes"=rownames(alt_table), "SNPs"=atac_snp, "w1"=w1,"w2"=w2)
  }
    rds_list[[paste0("chr",as.character(chrr))]]=result
    saveRDS(result,paste0(rds_path,"/EMresults/chr",as.character(chrr),".rds"))

    pdf(paste0(plot_path,"/EMresults/EMresult_chr", as.character(chrr),'.pdf' ))
    par(mfrow=c(2,1))
    hist(result$I_hat,100, xlim=c(0,1), main=paste0("Histogram of I_hat chr", as.character(chrr)))
    hist(result$theta_hat,100, xlim=c(0,1), main=paste0("Histogram of theta_hat chr", as.character(chrr)))
    dev.off()
    cat(paste0(chrr," "))
    }else{
      cat(paste0("skip ",chrr, " "))
      rds_list[[paste0("chr",as.character(chrr))]]=readRDS(paste0(rds_path,'/chr',chrr,'.rds'))
      next
    }
  }
  cat("\n")

  Obj_filtered$rds_list=rds_list

  message("Finsh iterations.")
  cat("Results and plots for each region have been saved in the rds and the plot directory.\n")
  cat("\"rds_list\" (EM result for each region) was added to the Obj_filtered object.")
  return(Obj_filtered)
}
