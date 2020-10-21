#' Genotype each cell based on coverage and major haplotype proportion (theta_hat) in a region.
#' 
#' possibleGenotypes, priorProb come from referece data.
#' 
#' @param X: A ncell by 2 dataframe. Column 1: normalized coverage (rho_hat); Column 2: theta_hat
#' @param snpCoverage: A numeric vector. Total read counts covering SNPs for each cell in a region.
#' @param possibleGenotypes: A numeric vector indicating all possible genotypes in a region.
#' @param priorProb: A numeric vector indicating prior prior probability of each possible genotype. The order should match the one for "possibleGenotypes".
#' @param sigma.rho: Numeric. Standard deviation of the rho_i values under normal distribution. 
#' 
#' @keywords internal
#' @export
genotype_ref<-function(X, snpCoverage, possibleGenotypes, priorProb,sigma.rho=0.25){
  rhohat=X[,1]
  thetahat=X[,2]
  canonicalPoints=cbind(c(1:27),
                        c(0.5,0.5, 1,1,1,1.5,1.5,1.5,1.5,2,2,2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3),
                        c(0,1,0,0.5,1,0,1/3,2/3,1,0,1/4,2/4,3/4,1,0,1/5,2/5,3/5,4/5,1,0,1/6,2/6,3/6,4/6,5/6,1))
  
  
  ncells=length(rhohat)
  ngenotypes=length(possibleGenotypes)
  genotypeProb = matrix(nrow=ncells, ncol=length(possibleGenotypes))
  genotypes=rep(NA, ncells)
  genotypeConfidence=rep(NA, ncells)
  colnames(genotypeProb)=paste(possibleGenotypes)
  for(i in 1:ncells){
    if(is.na(rhohat[i]) || is.na(thetahat[i])){
      genotypeProb[i,]=0
    } else {
      numerator=rep(0, ngenotypes)
      for(j in 1:ngenotypes){
        thetaj=canonicalPoints[possibleGenotypes[j], 3]
        if(thetaj==0){
          thetaj=0.05
        } 
        if(thetaj==1){
          thetaj=0.95
        }
        rhoj = canonicalPoints[possibleGenotypes[j], 2]
        numerator[j] = priorProb[j]*dnorm(rhohat[i], mean=rhoj, sd=sigma.rho)*dnorm(thetahat[i], mean=thetaj, sd=sqrt(thetaj*(1-thetaj)/snpCoverage[i]))
      }
      if(sum(numerator)==0){
        genotypeProb[i,]=0
        genotypes[i]=NA
      } else {
        genotypeProb[i,] = numerator/sum(numerator)
        genotypes[i]=possibleGenotypes[which.max(genotypeProb[i,])]
        genotypeConfidence[i]=max(genotypeProb[i,])
      }
    }
  }
  
  return(list(genotypes=genotypes, genotypeProb=genotypeProb, genotypeConfidence=genotypeConfidence))
}
