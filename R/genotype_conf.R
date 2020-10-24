#' Compute confidence scores based on posterior probability for each cell in a region.
#'
#' @param X: A ncell by 2 dataframe. Column 1: normalized coverage (rho_hat); Column 2: theta_hat
#' @param gt: A vector of lenth ncell. The numbers represent cell-level allele-specific copy number states.
#'
#' @return A lineage tree plot constructed using cell-level haplotype profiles across all regions.
#'
#' @import ggplot2
#' @import pheatmap
#' @import cluster
#' @export
genotype_conf=function(X=NULL, gt=NULL){
  
  canonicalPoints=cbind(c(1:27),
                        c(0.5,0.5, 1,1,1,1.5,1.5,1.5,1.5,2,2,2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3),
                        c(0,1,0,0.5,1,0,1/3,2/3,1,0,1/4,2/4,3/4,1,0,1/5,2/5,3/5,4/5,1,0,1/6,2/6,3/6,4/6,5/6,1))
  
  ncells=nrow(X)
  #nregions=ncol(Obj_filtered$genotypes)
  max.sdrho=rep(0.07, nrow(canonicalPoints))
  max.sdtheta=rep(0.06, nrow(canonicalPoints))
  max.sdrho[4] = 0.15
  max.sdtheta[4] = 0.1
  
  
  
  #posteriorConfidence=matrix(nrow=ncells, ncol=nregions, data=NA)
  posteriorConfidence=rep(0,ncells )
  
  #for(regionid in 1:nregions){
    #rhohat=Obj_filtered$genotype_values[,2*(regionid-1)+1]
    #thetahat=Obj_filtered$genotype_values[,2*(regionid-1)+2]
    #genotype=Obj_filtered$genotypes[,regionid]
  
  rhohat=X[,1]
  thetahat=X[,2]
  genotype=gt
    
    possible.genotypes=unique(genotype)
    cluster.props=rep(0, length(possible.genotypes))
    for(i in 1:length(possible.genotypes)) cluster.props[i]=sum(genotype==possible.genotypes[i])
    cluster.props=cluster.props/sum(cluster.props)
    
    mu.rho=rep(0, length(possible.genotypes))
    mu.theta=rep(0, length(possible.genotypes))
    sd.rho=rep(0, length(possible.genotypes))
    sd.theta=rep(0, length(possible.genotypes))
    for(i in 1:length(possible.genotypes)){
      cellids = which(genotype==possible.genotypes[i])
      mu.rho[i]=canonicalPoints[possible.genotypes[i],2]
      mu.theta[i]=canonicalPoints[possible.genotypes[i],3]
      sd.rho[i]=sd(rhohat[cellids])
      sd.theta[i]=sd(thetahat[cellids])
      if(is.na(sd.rho[i])) sd.rho[i] = max.sdrho[possible.genotypes[i]]
      if(is.na(sd.theta[i])) sd.theta[i] = max.sdtheta[possible.genotypes[i]]
      
    }
    sd.rho=pmin(sd.rho, max.sdrho[possible.genotypes])
    sd.theta=pmin(sd.theta,max.sdtheta[possible.genotypes])
    
    mm=match(genotype, possible.genotypes)
    numerator=rep(NA,ncells)
    denominator=rep(NA,ncells)
    for(i in 1:ncells){
      numerator[i]=cluster.props[mm[i]]*dnorm((rhohat[i]-mu.rho[mm[i]])/sd.rho[mm[i]])*dnorm((thetahat[i]-mu.theta[mm[i]])/sd.theta[mm[i]])
      denominator[i]=0
      for(j in 1:length(possible.genotypes)){
        denominator[i]=denominator[i]+cluster.props[j]*dnorm((rhohat[i]-mu.rho[j])/sd.rho[j])*dnorm((thetahat[i]-mu.theta[j])/sd.theta[j])
      }
      posteriorConfidence[i]=numerator[i]/denominator[i]  
    }
  #}
  
   return(posteriorConfidence)
  
}