#' Filter object based on cell number for each SNP, SNP number for each cell, SNP variant allele frequency, and exclude the centromere and telomere regions.
#'
#' @param Obj An Alleloscope object.
#' @param cell_filter An integer of minimum cell number for SNP selection.
#' @param SNP_filter An integer of minimum SNP number for cell selection.
#' @param min_vaf A numerical value in the range (0,1) of minimum SNP variant allele frequency in the pseudo bulk for SNP selection.
#' @param max_vaf A numerical value in the range (0,1) of mzsimum SNP variant allele frequency in the pseudo bulk for SNP selection.
#' @param centro A Matrix/ data.frame of centromere information.
#' @param telo A Matrix/ data.frame of telomere information.
#' @snp_ind A numeric vector indexing the SNPs to be included.
#' @param plot_stat Logical (TRUE/FALSE). Whether or not to plot the summary statistics.
#' @param plot_vaf Logical (TRUE/FALSE). Whether or not to plot the variant allele frequency for the pseudo bulk for all the chromosomes.
#'
#' @return A Alleloscope object after the filtering.
#'
#' @export
Matrix_filter=function(Obj=NULL, cell_filter=5, SNP_filter=10 ,min_vaf=0, max_vaf=1, centro=NULL, telo=NULL,snp_ind=NULL, plot_stat=TRUE, plot_vaf=TRUE){

  ## setting path and name
  alt_all=Obj$alt_all
  total_all=Obj$total_all
  var_list=Obj$var_all
  size=Obj$size
  assay=Obj$assay
  samplename=Obj$samplename
  plot_path=paste0(Obj$dir_path,'/plots/')
  
  af=Matrix::rowSums(alt_all)/Matrix::rowSums(total_all)
  af[is.na(af)]=0

  ## plot stat
  if(plot_stat==TRUE){
    pdf(paste0(plot_path,"statistics.pdf" ))
    par(mfrow=c(3,1))
    hist(Matrix::colSums(total_all), main=paste0(as.character(dim(total_all)[2])," cells across ",  as.character(dim(total_all)[1]), " SNPs"), xlab="coverage of individul cells", ylab='frequency', breaks = 100)
    hist(Matrix::rowSums(total_all), main=paste0(as.character(dim(total_all)[1])," SNPs across ",  as.character(dim(total_all)[2]), " cells"), xlab="coverage of individul SNPs", xlim=c(0,100), ylab='frequency', breaks = 10000)
    hist(af, 1000, main="Histogram of VAF values")
    dev.off()
  }
  
  ##
  if(!is.null(snp_ind)){
    if(!is.numeric(snp_ind)){
      message("Please provide numeric snp_ind for the set of SNPs you want to include.")
    }else{
    alt_all=alt_all[snp_ind,]
    total_all=total_all[snp_ind,]
    var_list=var_list[snp_ind,]}
  }
  
  
  cc=Matrix::colSums(alt_all)
  cc_ind=which(cc>cell_filter)

  message(paste0(length(cc_ind), " cells after filtering."))

  total_all=total_all[,(cc_ind) ]
  alt_all=alt_all[,(cc_ind)]
  af=Matrix::rowSums(alt_all)/Matrix::rowSums(total_all)
  af[is.na(af)]=0


  rr=Matrix::rowSums(total_all)
  if(assay=="scDNAseq"){
  rr_ind=which(rr>SNP_filter & rr<= ((median(rr[which(rr!=0)])+3*mad(rr[which(rr!=0)]))) & af<=max_vaf & af>min_vaf) ###
  }else{
  rr_ind=which(rr>SNP_filter & af<=max_vaf & af>min_vaf) }
  
  
  message(paste0(length(rr_ind), " SNPs after filtering."))
  cat("(Recommend more than 1,000,000 SNPs for all chromosomes)\n")

  total_all=total_all[(rr_ind), ]
  alt_all=alt_all[(rr_ind),]

  ## calcuate new VAF for the filtered matrices
  af=Matrix::rowSums(alt_all)/Matrix::rowSums(total_all)
  af[is.na(af)]=0

  # var_list
  var_list=var_list[rr_ind,]
  if(sapply(strsplit(as.character(var_list[2,1]),'hr'),'[',1)=='c'){
    var_list[,1]=sapply(strsplit(var_list[,1],'hr'),'[',2) ## if the varlist has 'chr'
  }
  ## check the order of the SNPs
  oo=order(as.numeric(as.character(var_list[,1])), as.numeric(as.character(var_list[,2])))
  var_list=var_list[oo,]
  af=af[oo]
  total_all=total_all[oo,]
  alt_all=alt_all[oo,]

  ## get rid of telomere and centromere

  if(!(is.null(centro)) | !(is.null(telo))){
  if(!(is.null(centro)) & (is.null(telo))){
    message("Filter SNPs in the centromere regions.")
    centro$V2=sapply(strsplit(centro$V2,'hr'),'[',2)
    centro_tele=centro[,2:4]
  }else if((is.null(centro)) & !(is.null(telo))){
    message("Filter SNPs in the telomere regions.")
    telo$V2=sapply(strsplit(telo$V2,'hr'),'[',2)
    centro_tele=telo[,2:4]
  }else if(!(is.null(centro)) & !(is.null(telo))){
    message("Filter SNPs in the centromere and telomere regions.")
    telo$V2=sapply(strsplit(telo$V2,'hr'),'[',2)
    centro$V2=sapply(strsplit(centro$V2,'hr'),'[',2)
    centro_tele=rbind(centro[,2:4], telo[,2:4])
  }


  ## remove centromeres and telomeres
  if(assay=="scDNAseq"){
    rm_list=c()
    for(ii in 1:22){
      centro_tele_sub=centro_tele[which(centro_tele[,1]==as.character(ii)),]
      for(jj in 1:dim(centro_tele_sub)[1]){
        rm_ind=which(var_list[,1]==as.character(ii) & (as.numeric(var_list[,2])>= as.numeric(centro_tele_sub[jj,2]) & as.numeric(var_list[,2])<= as.numeric(centro_tele_sub[jj,3])) )
        rm_list=c(rm_list, rm_ind)}
    }
    rm_list=unique(rm_list)
    if(length(rm_list)>0){
    var_list=var_list[-rm_list,]
    af=af[-rm_list]
    total_all=total_all[-rm_list,]
    alt_all=alt_all[-rm_list,]}
    message(paste0(as.character(length(rr_ind)-length(rm_list)), " SNPs after centromeres and telomeres filtering."))
  }
  }

  var_str=paste0(as.character(var_list[,1]),":", as.character(var_list[,2]),"_", as.character(var_list[,4]),"_", as.character(var_list[,5]))
  var_chr=as.numeric(sapply(strsplit(var_str,':'),'[',1))
  var_pos=as.numeric(sapply(strsplit(sapply(strsplit(var_str,':'),'[',2), "_"),'[',1))
  chr_nvar=rep(0, length(size))
  chr_nvar=as.numeric(table(var_chr))



  ### for segmentation
  rownames(total_all)=var_str
  rownames(alt_all)=var_str
  #ref_all=total_all-alt_all

  #############

  size_cs=c(0,cumsum(size))
  size_cs_rep=rep(size_cs[1:(length(size_cs)-1)], chr_nvar)
  loci_xaxis=var_pos+size_cs_rep
  chrx=size/2+size_cs[1:length(size)]


  ## plot total chromosome
  if(plot_vaf==TRUE){
    if(length(af)<50000){
      ind=1:length(af)
    }else{
      ind=sort(sample(1:length(af),50000, replace = F))
    }
    pdf(paste0(plot_path,"AF_",assay,"_c",as.character(cell_filter),"_r", as.character(SNP_filter),".pdf"), width = 24, height = 12)
    par(mfrow=c(1,1))
    plot(loci_xaxis[ind], af[ind], xlab="", ylab='AF', main=paste0('Allele frequency ',samplename,'_', assay), xaxt="n", yaxt="n", xlim = c(0, size_cs[length(size_cs)]), ylim=c(0,1),pch=20, cex=1 , cex.lab=2, cex.main=3)
    axis(side=1, at=chrx, labels = names(size), cex.axis=2)
    axis(side=2, at=c(0, 0.5,1), labels = c(0, 0.5, 1), cex.axis=2)
    abline(v=size_cs, col='red', lty=1, lwd=4)
    dev.off()
  }

  filter=list("alt_all"=alt_all, "total_all"=total_all,"var_all"=var_list,
              'cell_filter'=cell_filter, 'SNP_filter'=SNP_filter ,'min_vaf'=min_vaf, 'max_vaf'=max_vaf,"barcodes"=Obj$barcodes,
              "size"=Obj$size, 'samplename'=Obj$samplename, 'dir_path'=Obj$dir_path,
              'genome_assembly'=Obj$genome_assembly, 'assay'=Obj$assay)



  message("Object successfully filterd!")
  cat(paste0("Plots for statistics have been saved in the path:", dir_path,"plots/",'\n'))
  return(filter)
}
