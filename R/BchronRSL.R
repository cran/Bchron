BchronRSL <-
function(Bchrondata,RSLmean,RSLsd,degree=1,iter=10000,burnin=2000,thin=8,reportevery=1000,BP=FALSE) {
  #Bchrondata=mycore;RSLmean=IcelandRSL[,2];RSLsd=IcelandRSL[,3];degree=1;iter=10000;burnin=2000;thin=8;reportevery=1000;BP=FALSE
  
  if(degree>5) stop('Degree not supported')
  pow.names = c('mean','linear','quadratic','cubic','quartic','quintic')  
  remaining=(iter-burnin)/thin
  beta.store = matrix(NA,ncol=degree+1,nrow=remaining)
  y = matrix(RSLmean,ncol=1)
  Q = diag(1/((RSLsd)^2))
  N = nrow(y)
  chrons = read.table(Bchrondata$chrons)
  whichrows = sample(1:nrow(chrons),iter,replace=TRUE)
  degmat = matrix(rep(0:(degree),ncol(chrons)*N),nrow=N,ncol=degree+1,byrow=TRUE)
  if(BP) const = mean(as.matrix(chrons))
  if(!BP) const = mean(1950-1000*as.matrix(chrons))
  for(i in 1:iter) {
    if(i%%reportevery==0) cat('\r',round(100*i/iter,1),'%')
    if(i%%20==0 | i==1) {
      if(BP) currchron = chrons[whichrows[i],] - const
      if(!BP) currchron = 1950-1000*chrons[whichrows[i],] - const
      X = matrix(rep(as.numeric(currchron),degree+1),ncol=degree+1)
      X = X^degmat
    }
    # Sample a beta
    if(i%%thin==0 & i>burnin) {
      beta.store[(i-burnin)/thin,] = matrix(rmvnorm(1,solve(t(X)%*%Q%*%X,t(X)%*%Q%*%y),solve(t(X)%*%Q%*%X)),ncol=1)
    }
  }
  
  cat('\r')
  cat('95% posterior intervals...\n')
  cat('Power Lower Upper \n')
  for(j in 1:(degree+1)) {
    cat(pow.names[j],quantile(beta.store[,j],probs=c(0.025,0.975)),'\n')
  }
  
  return(list(samples=beta.store,degree=degree,BP=BP,RSLmean=RSLmean,RSLsd=RSLsd,chrons=chrons,const=const))
  
}
