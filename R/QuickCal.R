Quickcal <- function(date,sd,pathtocalcurve=NULL,prob=c(68,95,99),get.output=FALSE) {
  
  calcurve = as.matrix(read.table(pathtocalcurve,sep=','))
  calbp = calcurve[,1]
  c14bp = calcurve[,2]
  calsd = calcurve[,3]
  
  # Get rid of dates outside the range of the calibration curve
  if(date>max(calbp) | date<min(calbp)) stop("Radiocarbon date outside of calibration curve range")
  
  # Create an age grid and get mean and variance of calibration curve 
  agegrid = seq(min(calbp),max(calbp),by=1)
  mu = approx(calbp,c14bp,xout=agegrid)$y
  tau = sd^2 + approx(calbp,calsd,xout=agegrid)$y
  
  # Calculate density
  dens = dnorm(date,mean=mu,sd=sqrt(tau))
  dens = dens/sum(dens)
  
  # Get hdr and print
  out = list(x=agegrid[dens>0],y=dens[dens>0])
  cat('Highest density regions for date',date,'+/-',sd,'\n')
  print(round(hdr(den=out,prob=prob)$hdr,1))
  
  if(get.output) return(list(age=out$x,density=out$y))
  
}