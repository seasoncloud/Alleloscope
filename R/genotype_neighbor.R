#' Genotype each cell based on coverage and major haplotype proportion (theta_hat) in a region.
#'
#' @param X: A ncell by 2 dataframe. Column 1: normalized coverage (rho_hat); Column 2: theta_hat
#' @noRd
#' @keywords internal
#' @export
genotype_neighbor=function(X, cluster_name=F){
  mu=matrix(c(c(0.5,0), c(0.5,1),c(1,0),c(1,0.5),c(1,1),c(1.5,0),c(1.5,0.33),c(1.5, 0.66),c(1.5, 1),c(2, 0), c(2, 0.25), c(2,0.5),c(2,0.75),c(2,1),c(2.5, 0),c(2.5,0.2),c(2.5,0.4), c(2.5, 0.6), c(2.5, 0.8), c(2.5, 1), c(3,0),c(3,1/6),c(3,2/6),c(3,3/6),c(3,4/6),c(3,5/6),c(3,1))
            , byrow=T, ncol=2)

  center_dict=1:27
  names(center_dict)=paste0("rho",mu[,1],"_theta", round(mu[,2],2))

  d_matrix=matrix(nrow=dim(X)[1], ncol=dim(mu)[1])
  for(ii in 1:dim(mu)[1]){
    dd=apply(X,1,function(ss) sqrt(sum((ss-mu[ii,])^2)))
    d_matrix[,ii]=dd}

  cluster=apply(d_matrix,1, function(ss) which(ss==min(ss))[1])
  sel_cluster=as.numeric(names(table(cluster))[which(table(cluster)>0)])

  ## update
  mu=mu[sel_cluster,, drop=F]
  d_matrix=matrix(nrow=dim(X)[1], ncol=dim(mu)[1])
  for(ii in 1:dim(mu)[1]){
    dd=apply(X,1,function(ss) sqrt(sum((ss-mu[ii,])^2)))
    d_matrix[,ii]=dd}
  cluster=apply(d_matrix,1, function(ss) which(ss==min(ss))[1])
  names(cluster)=paste0("rho",mu[cluster,1],"_theta", mu[cluster,2])
  if(cluster_name==FALSE){
  cluster_ind=center_dict[names(cluster)]
  }else{
  cluster_ind=names(cluster)
  }
  return(cluster_ind)
}
