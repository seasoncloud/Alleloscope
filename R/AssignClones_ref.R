#' Using marker regions to assign each cell into c reference subclones
#' 
#' rhohats, thetahats, snpCoverages are n by m matrices for n cell and m marker regions.
#' The genotypes (rho, theta) are: 1.(0.5,0); 2.(0.5,1); 3.(1,0); 4.(1,0.5); 5.(1,1); 6.(1.5,0); 7.(1.5,1/3); 8.(1.5,2/3); 9.(1.5,1); 10.(2,0); 11.(2,1/4); 12.(2,2/4); 13.(2,3/4); 14.(2,4/4)
#'                                15.(2.5,0); 16.(2.5,1/5); 17.(2.5,2/5); 18.(2.5,3/5); 19.(2.5,4/5); 20.(2.5,1); 21.(3,0); 22.(3,1/6); 23.(3,2/6); 24.(3,3/6); 25.(3,4/6); 26.(3,5/6); 27.(3,6/6)
#' @param Obj_filtered An Alleloscope object with a n cell by (m region * 4) genotype_values matrix.
#' Every 4 columns in the genotype_table matrix are (rho_hat, theta_hat, h1, h2) of each region.
#' h1, h2 are the coverage across all SNPs located on the major haplotype (h1) and the minor haplotype (h2) in a region for each cell
#' @param priorCloneProbs: A numeric vector indicating prior prior probability of each subclone.
#' @param clone.genotypes: c by m matrix of numbers representing different genotypes for each clone and each maker region (known from scDNA-seq). 
#' @param sigma.rho: Numeric. Standard deviation of the rho_i values under normal distribution. 
#' 
#' @export
AssignClones_ref<-function(Obj_filtered=NULL, priorCloneProbs=NULL, clone.genotypes, sigma.rho=0.25){
  atacdat=Obj_filtered$genotype_values
  ncells=nrow(atacdat)
  #nregions=length(markerRegions)
  nregions = ncol(clone.genotypes)
  nclones=nrow(clone.genotypes)
  
  rhohats=matrix(nrow=ncells, ncol=nregions)
  thetahats=matrix(nrow=ncells, ncol=nregions)
  snpCoverages=matrix(nrow=ncells, ncol=nregions)
  for(j in 1:nregions){
    regionid=j
    rhohats[,j]=atacdat[,(regionid-1)*4+1]
    thetahats[,j]=atacdat[,(regionid-1)*4+2]
    h1=atacdat[,(regionid-1)*4+3]
    h2=atacdat[,(regionid-1)*4+4]
    snpCoverages[,j]=h1+h2
  }
  
  if(is.null(priorCloneProbs)){
    priorCloneProbs=rep(1/6,6)
  }
  canonicalPoints=cbind(c(1:27),
                        c(0.5,0.5, 1,1,1,1.5,1.5,1.5,1.5,2,2,2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3),
                        c(0,1,0,0.5,1,0,1/3,2/3,1,0,1/4,2/4,3/4,1,0,1/5,2/5,3/5,4/5,1,0,1/6,2/6,3/6,4/6,5/6,1))
  
  #nregions = ncol(clone.genotypes)
  #nclones=nrow(clone.genotypes)
  #ncells=nrow(rhohats)
  theta.kjs = matrix(nrow=nclones, ncol=nregions,canonicalPoints[clone.genotypes,3], byrow=FALSE)
  rho.kjs= matrix(nrow=nclones, ncol=nregions,data=canonicalPoints[clone.genotypes,2], byrow=FALSE)
  theta.kjs=pmax(theta.kjs, 0.05)
  theta.kjs=pmin(theta.kjs, 0.95)
  
  cloneLogLikelihoods=matrix(nrow=ncells, ncol=nclones)
  probClone=matrix(nrow=ncells, ncol=nclones)
  cloneAssign=rep(NA,ncells)
  cloneConf=rep(NA, ncells)
  for(i in 1:ncells){
    for(k in 1:nclones){
      notNA=which(!is.na(rhohats[i,]) & !is.na(thetahats[i,]) & !is.na(snpCoverages[i,]))
      sdtheta=sqrt(theta.kjs[k,notNA]*(1-theta.kjs[k,notNA])/snpCoverages[i,notNA])
      cloneLogLikelihoods[i,k]=sum(-log(sqrt(2*pi)*sigma.rho)-(rhohats[i,notNA]-rho.kjs[k,notNA])^2/(2*sigma.rho^2))+
        sum(-log(sqrt(2*pi)*sdtheta)-(thetahats[i,notNA]-theta.kjs[k,notNA])^2/(2*sdtheta^2))
    }
    minl = min(cloneLogLikelihoods[i,])
    for(k in 1:nclones){
      probClone[i,k]=exp(cloneLogLikelihoods[i,k]-minl)*priorCloneProbs[k]
    }
    probClone[i,]=probClone[i,]/sum(probClone[i,])
    sel=which(!is.na(probClone[i,]))
    if(length(sel)>0){
      cloneAssign[i]=which.max(probClone[i,])
      cloneConf[i]=max(probClone[i,], na.rm=TRUE)
    } else {
      cloneAssign[i]=NA
      cloneConf[i]=NA
      
    }
  }
    cloneAssign=list(probClone=probClone, cloneAssign=cloneAssign, cloneConf=cloneConf)
    Obj_filtered$cloneAssign=cloneAssign
  return(Obj_filtered)
}