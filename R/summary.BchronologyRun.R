# Summary method for BchronologyRun
summary.BchronologyRun = function(object,probs=c(0.025,0.1,0.5,0.9,0.975), ..., digits = max(3, getOption("digits")-3)) {
  # Give a list of quantiles for each depth and Show outlier probabilities
  chronSummary = data.frame(cbind(object$predictPositions,round(t(apply(object$thetaPredict,2,'quantile',probs)),3)))
  colnames(chronSummary) = c('Depth',paste(probs*100,'%',sep=''))
  cat('Quantiles of predicted ages by depth: \n')
  print(chronSummary,row.names=FALSE)
  cat('\n')
  cat('Posterior outlier probability by date: \n')
  outprob = data.frame(names(object$calAges),colMeans(object$phi))
  colnames(outprob) = c('Date','OutlierProb')
  print(outprob,row.names=FALSE)
  cat('\n')
  cat('Convergence check (watch for too many small p-values): \n')
  pars = cbind(object$theta,object$phi,object$mu,object$psi)
  n = ncol(object$theta)
  colnames(pars) = c(names(object$calAges),paste('Outlier',1:n),'RateMean','RateVar')
  geweke = coda::geweke.diag(pars)[[1]]
  geweke[is.nan(geweke)] = 0
  pvals = data.frame(sort(round(c(pnorm(geweke[geweke<0]),1-pnorm(geweke[geweke>0])),5)))
  colnames(pvals) = 'p-value'
  print(pvals)  
}