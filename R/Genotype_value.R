#' Normalize coverage using identified/ specified normal cells and one normal region and generate a table with (rho_hat, theta_hat) of each cell for all regions.
#'
#' rho_hat: Relative coverage change for each cell in a region
#' theta_hat: Major haplotype proportion fir each cell in a region
#'
#' @param Obj_filtered An Alleloscope object with theta_hat info in the rds_list and identified/ specified normal cells and a normal region
#' @param type Specify whether the sample is a "tumor" or "cellline". If "type" is a "cellline", param "ref_counts" needs to be specified for normal sample.
#' @param raw_counts (required) A large binned coverage matrix (m1 bin by n1 cell) with values being read counts for all chromosomal regions of tumor sample.
#' @param ref_counts (required only when type = "cellline") A binned coverage matrix (m2 bin by n2 cell) with values being read counts for all chromosomal regions of normal sample. n2 can be 1 for bulk sample.
#' @param cov_adj An integer for coverage adjustment for tumor cells. 
#' @param ref_gtv A reference "genotype_values" (from scDNA-seq) to help with rho_i estimation.
#' @param mincell An integer to filter out regions with minimum number of cells.
#' @param cell_filter Logical (TRUE/ FALSE). Whether or not to exclude cells with rho_hat>0.99 or <0.01 for each region. 
#' @param cell_filter Logical (TRUE/ FALSE). Whether or not to exclude low quality cells in the output matrix. 
#' @param refr Logical (TRUE/ FALSE). Whether or not to use diplid region for normalization (otherwise, cell size is used).
#' 
#' @return (rho_hat, theta_hat) of each cell for all region in the "genotype_values".
#' Every 2 columns in the genotype_table are (rho_hat, theta_hat) of each region. Each row is a cell.
#'
#' @export
Genotype_value=function(Obj_filtered=NULL, type="tumor", raw_counts=NULL, ref_counts=NULL, cov_adj=1, ref_gtv=NULL, mincell=NULL, qt_filter=TRUE,  cell_filter=TRUE, refr=TRUE){
  samplename=Obj_filtered$samplename
  dir_path=Obj_filtered$dir_path
  assay=Obj_filtered$assay
  barcode_normal=Obj_filtered$select_normal$barcode_normal
  ref=Obj_filtered$ref
  refn=sapply((strsplit(ref,":")),'[',1)
  #cell_info=Obj_filtered$cell_info
  seg_table_filtered=Obj_filtered$seg_table_filtered
  #seg_table=Obj_filtered$seg_table
  rds_list=Obj_filtered$rds_list
  cell_barcodes=Obj_filtered$barcodes

  ncell=length(Obj_filtered$barcodes)
  theta_N_nr_nc=list()

  
  ##ref
  if(!is.null(ref_gtv)){
  dna_gt=ref_gtv
  dna_gt=dna_gt[,which(sapply(strsplit(colnames(dna_gt),"_"),'[',1)=='rho')]
  }
  
  if(is.null(mincell)){
    mincell=ncol(Obj_filtered$total_all)*0.9
  }
  
  ## raw/ref count matrix info
  raw_chr=sapply(strsplit(rownames(raw_counts),'-'),'[',1)
  raw_start=as.numeric(sapply(strsplit(rownames(raw_counts),'-'),'[',2))
  raw_end=as.numeric(sapply(strsplit(rownames(raw_counts),'-'),'[',3))
  
  if(is.null(ref_gtv) & type=='cellline'){
    if(refr==TRUE){
  ref_chr=sapply(strsplit(rownames(ref_counts),'-'),'[',1)
  ref_start=as.numeric(sapply(strsplit(rownames(ref_counts),'-'),'[',2))
  ref_end=as.numeric(sapply(strsplit(rownames(ref_counts),'-'),'[',3))}
  }
  
  N0_all=colSums(raw_counts[,match(result$barcodes, cell_barcodes)]) ## for p and q #for cytoarm
  
  message("Start estimating cell specific (rho_hat, theta_hat) for each region.")


  
  ####
  regions=gsub('chr','',names(Obj_filtered$rds_list))
  for(chrr in regions){
    chrrn=unlist(strsplit(chrr,':'))[1]

    result=rds_list[[paste0('chr',chrr)]]
    theta_hat=result$theta_hat
    names(theta_hat)=result$barcodes

    raw_counts_chr=raw_counts[which(raw_chr %in% paste0('chr', as.character(chrrn))),]
    raw_chr_sub=raw_chr[which(raw_chr %in% paste0('chr', as.character(chrrn)))]
    raw_start_sub=raw_start[which(raw_chr %in% paste0('chr', as.character(chrrn)))]
    raw_end_sub=raw_end[which(raw_chr %in% paste0('chr', as.character(chrrn)))]
    
    
    
    
    ## subsetting region
    query=GenomicRanges::GRanges(paste0('chr',chrrn),IRanges::IRanges(as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == chrr),2]), as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == chrr),3])))
    subject=GenomicRanges::GRanges(raw_chr_sub, IRanges::IRanges(as.numeric(raw_start_sub),as.numeric(raw_end_sub))) ## cytoband 1-based start and 1-based end
    ov=findOverlaps(query, subject)
    ov=as.matrix(ov)
    
    bin_start=min(ov[,2])
    bin_end=max(ov[,2])
    raw_counts_chr=raw_counts_chr[bin_start:bin_end,match(result$barcodes, cell_barcodes)] ## for p and q #for cytoarm
    Nr=colSums(raw_counts_chr)
    names(Nr)=result$barcodes
    

    # subsetting normal
    if(refr==TRUE){
      raw_counts_ref=raw_counts[which(raw_chr %in% paste0('chr', as.character(refn))),]
      raw_ref_chr_sub=raw_chr[which(raw_chr %in% paste0('chr', as.character(refn)))]
      raw_ref_start_sub=raw_start[which(raw_chr %in% paste0('chr', as.character(refn)))]
      raw_ref_end_sub=raw_end[which(raw_chr %in% paste0('chr', as.character(refn)))]
    
    query=GenomicRanges::GRanges(paste0('chr',refn),IRanges::IRanges(as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == ref),2]), as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == ref),3])))
    subject=GenomicRanges::GRanges(raw_ref_chr_sub, IRanges::IRanges(as.numeric(raw_ref_start_sub),as.numeric(raw_ref_end_sub))) ## cytoband 1-based start and 1-based end
    ov=findOverlaps(query, subject)
    ov=as.matrix(ov)
    
    bin_start=min(ov[,2])
    bin_end=max(ov[,2])
    raw_counts_ref=raw_counts_ref[bin_start:bin_end,match(result$barcodes, cell_barcodes)] ## for p and q #for cytoarm
    N0=colSums(raw_counts_ref)
    }else{
      #raw_counts_ref=raw_counts[,match(result$barcodes, cell_barcodes)] ## for p and q #for cytoarm
      N0=N0_all
      ref="cell_size"
      Obj_filtered$ref=ref
    }

    ## non_noisy
    barcodes_non_noisy=cell_barcodes#cell_info$barcode[which(cell_info$is_noisy==0)]
    barcodes_non_noisy=intersect(barcodes_non_noisy, result$barcodes)


    Nr=Nr[match(barcodes_non_noisy, result$barcodes)]
    N0=N0[match(barcodes_non_noisy, result$barcodes)]
    Ni=Nr/N0
    names(Ni)=barcodes_non_noisy
    
    if(is.null(ref_gtv)){
      if(length(result$barcodes)<mincell){
        cat(paste0("Exclude ",chrr," region:<",mincell," cells\n"))
        next
      }
    if(type=='tumor'){
      ref_ncell=length(barcode_normal)
      Nrref=Nr[which(names(Nr) %in% barcode_normal)]
      N0ref=N0[which(names(Nr) %in% barcode_normal)]
      Ni_ref=median(Nrref/N0ref)

    }else if(type=='cellline'){
      #normal sample from other normal dataset
        ref_counts_chr=ref_counts[which(ref_chr %in% paste0('chr', as.character(chrrn))),]
        ref_counts_ref=ref_counts[which(ref_chr %in% paste0('chr', as.character(refn))),]
        
        ref_chr_sub=ref_chr[which(ref_chr %in% paste0('chr', as.character(chrrn)))]
        ref_start_sub=ref_start[which(ref_chr %in% paste0('chr', as.character(chrrn)))]
        ref_end_sub=ref_end[which(ref_chr %in% paste0('chr', as.character(chrrn)))]
        
        ref_ref_chr_sub=ref_chr[which(ref_chr %in% paste0('chr', as.character(refn)))]
        ref_ref_start_sub=ref_start[which(ref_chr %in% paste0('chr', as.character(refn)))]
        ref_ref_end_sub=ref_end[which(ref_chr %in% paste0('chr', as.character(refn)))]
        
        # ref chr
        query=GenomicRanges::GRanges(paste0('chr',chrrn),IRanges::IRanges(as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == chrr),2]), as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == chrr),3])))
        subject=GenomicRanges::GRanges(ref_chr_sub, IRanges::IRanges(as.numeric(ref_start_sub),as.numeric(ref_end_sub))) ## cytoband 1-based start and 1-based end
        ov=findOverlaps(query, subject)
        ov=as.matrix(ov)
        
        bin_start=min(ov[,2])
        bin_end=max(ov[,2])
        ref_counts_chr=ref_counts_chr[bin_start:bin_end,] ## for p and q #for cytoarm
        
        ## ref ref
        query=GenomicRanges::GRanges(paste0('chr',refn),IRanges::IRanges(as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == ref),2]), as.numeric(seg_table_filtered[which(seg_table_filtered$chrr == ref),3])))
        subject=GenomicRanges::GRanges(ref_ref_chr_sub, IRanges::IRanges(as.numeric(ref_ref_start_sub),as.numeric(ref_ref_end_sub))) ## cytoband 1-based start and 1-based end
        ov=findOverlaps(query, subject)
        ov=as.matrix(ov)
        
        bin_start=min(ov[,2])
        bin_end=max(ov[,2])
        ref_counts_ref=ref_counts_ref[bin_start:bin_end,] ## for p and q #for cytoarm
        
      Nrref=colSums(ref_counts_chr)
      N0ref=colSums(ref_counts_ref)
      #
      Ni_ref=Nrref/N0ref
      Ni_ref[is.na(Ni_ref)]=0
      Ni_ref=median(Ni_ref)
    }
    
  if(type=='cellline'){
    Ni=Ni*cov_adj
  }else{
    Ni[which(! names(Nr) %in% barcode_normal)]=Ni[which(! names(Nr) %in% barcode_normal)]*cov_adj
    }
    
    Ni=Ni/Ni_ref
  if(qt_filter==TRUE){
    Niq=Ni[which(Ni<=quantile(Ni, 0.99) & Ni>=quantile(Ni, 0.01))] ##
  }else{
    Niq=Ni
    Niq=pmin(Niq, quantile(Ni, 0.99))
    Niq=pmax(Niq, quantile(Ni, 0.01))
  }
    barcodes_nn_q=names(Niq)
    
    }else{
      gt=as.numeric(dna_gt[,which(colnames(dna_gt)==paste0('rho_',chrr))])
      if(length(Ni)>(ncell/2)){
        md=median(gt)
        Ni=md/median(Ni, na.rm = TRUE)*Ni
      }else{
        md=quantile(gt,0.2)
        Ni=md/quantile(Ni,0.2, na.rm=TRUE)*Ni}
      Niq=Ni
      barcodes_nn_q=names(Niq)
      
      
    }

    theta_hat=theta_hat[match(barcodes_nn_q, names(theta_hat))]

    w1=result$w1[match(names(Ni), names(result$w1))]
    w2=result$w2[match(names(Ni), names(result$w2))]
    
    theta_N_nr_nc[[paste0("rho_",as.character(chrr))]]=Niq
    theta_N_nr_nc[[paste0("theta_",as.character(chrr))]]=theta_hat
    theta_N_nr_nc[[paste0("h1_",as.character(chrr))]]=w1
    theta_N_nr_nc[[paste0("h2_",as.character(chrr))]]=w2
    

    cat(paste0(chrr," "))
  }
  cat('\n')

  cell_list<-lapply(theta_N_nr_nc, function(x) {
    names(x)
  })

  if(!is.null(ref_gtv)){
    cell_intersect <- Reduce(union, cell_list)
  }else{
    if(cell_filter==TRUE){
      cell_intersect <- Reduce(intersect, cell_list)
    }else{
      cell_intersect <- Reduce(union, cell_list) 
    }
  }


  theta_hat_cbn <- sapply(theta_N_nr_nc,function(x){
    x[match(cell_intersect, names(x))]
  })

  rownames(theta_hat_cbn) <- cell_intersect

  saveRDS(theta_hat_cbn, paste0(Obj_filtered$dir_path, "/rds/genotype_values.rds"))

  Obj_filtered$genotype_values=theta_hat_cbn
  message("\"genotypes_values\" is added to the Obj_filtered object.")
  cat(paste0("Matrix for cell specific (theta, rho) for each region is stored as genotype_values.rds in the path:",dir_path,"\n"))




  return(Obj_filtered)
}
