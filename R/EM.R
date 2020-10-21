#' Iterative phasing and theta_hat estimation
#'
#' @param ref_table A SNP by cell read count matrix/ spare matrix for the reference alleles.
#' @param alt_table A SNP by cell read count matrix/ spare matrix for the alternative alleles.
#' @param max_iter An integer of maximum iteration number.
#' @param seed An integer of random seed number for EM initialization.
#'
#' @return A list of estimated indicators (I_hat) for each SNP and estimated major haplotype proportion (theta_hat) for each cell in one region. I_hat is the phasing result indicating whether reference allele is on the major haplotype for each SNP. Theta_hat represents the CNV states for each cell. A cell is considered as a CNV carrier if its theta_hat depart from 0.5.
#'
#' @export
EM=function(ref_table, alt_table, max_iter=max_iter, seed = 2020){
  alt_table=t(alt_table)
  ref_table=t(ref_table)
  tot_table=alt_table+ref_table
  mm=dim(ref_table)[2]  # mm snv
  nn=dim(ref_table)[1]  # nn cell


  var_tot=colSums(tot_table)
  var_alt=colSums(alt_table)
  var_vaf=var_alt/var_tot


  #k means clustering to get priors
  set.seed(seed)
  km=kmeans(x = var_vaf, centers = 3)
  km_label=rep(0.5, mm)
  oo=order(km$centers)
  km_label[which(km$cluster==oo[1])]=km$centers[oo[3]]
  km_label[which(km$cluster==oo[3])]=km$centers[oo[1]]

  ## EM

  ind0=km_label
  ind=ind0

  tol=0.001
  ll_old=-Inf


  for(ii in 1:max_iter){
    ind_table=matrix(rep(ind, nn), nrow = nn, byrow = T)
    # maximization step
    w1=rowSums((ref_table*ind_table)+(alt_table*(1-ind_table)))
    w2=rowSums((ref_table*(1-ind_table))+(alt_table*ind_table))

    theta=w1/(w1+w2)
    theta[is.na(theta)]=0.5

    ## estimation step
    theta_table=matrix(rep(theta, mm),nrow=nn, byrow=F)
    product= matrixStats::colProds(((theta_table)/(1-theta_table+1e-10))^(alt_table-ref_table))
    ind=1/(1+product)
    ll_new= sum(log(theta+1e-10)*w1 +log(1-theta+1e-10)*w2)


    if(abs(ll_new-ll_old)<tol){
      break}
    ll_old=ll_new
  }
  #dev.off()
  output=list(theta_hat=theta, I_hat=ind, iterations=ii, w1=w1, w2=w2)
  return(output)
}
