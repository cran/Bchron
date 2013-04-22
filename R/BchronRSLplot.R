BchronRSLplot <-
function(Bchrondata,RSLresults,xgrid=NULL,colours=colorRampPalette(c('white','black','black','black','black')),transp=0.2,ss=1000,ellipses=TRUE,blobs=FALSE) {
  #Bchrondata=mycore;RSLresults=RSLresults;xgrid=NULL;colours=colorRampPalette(c('white','black','black','black','black'));transp=0.2;ss=1000;ellipses=TRUE
  
  ranges = read.table(Bchrondata$rangesfile,header=TRUE)
  if(RSLresults$BP) {
    age.med = ranges[,3]
    age.low = ranges[,2]
    age.high = ranges[,4]    
  } else {
    age.med = 1950-1000*ranges[,3]
    age.low = 1950-1000*ranges[,2]
    age.high = 1950-1000*ranges[,4]
  }
  
  RSL.low = RSLresults$RSLmean-2*RSLresults$RSLsd
  RSL.high = RSLresults$RSLmean+2*RSLresults$RSLsd
  
  par(mar=c(4,4,3,1))
  if(RSLresults$BP) plot(age.med,RSLresults$RSL.mean,xlab='Age (k BP)',ylab='RSL (m)',type='n',las=1,xlim=range(c(age.low,age.high)),ylim=range(c(RSL.low,RSL.high)))
  if(!RSLresults$BP) plot(age.med,RSLresults$RSL.mean,xlab='Year AD',ylab='RSL (m)',type='n',las=1,xlim=range(c(age.low,age.high)),ylim=range(c(RSL.low,RSL.high)))
  title(Bchrondata$fullname)
  
  add.alpha <- function(COLOURS, MINALPHA){
    if(missing(MINALPHA)) stop("provide a value for alpha between 0 and 1")
    RGB <- col2rgb(COLOURS, alpha=TRUE)
    ALPHASEQ = seq(MINALPHA,1,length=ncol(RGB))
    RGB[4,] = round(RGB[4,]*ALPHASEQ)
    NEW.COLOURS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
    return(NEW.COLOURS)
  }
  pal <- colours
  COLOURS <- add.alpha(pal(20), transp)
  
  if(blobs) {
    for(i in 1:length(RSLresults$RSLmean)) {
      # Create a bivariate density plot and then create some kind of contour
      if(!RSLresults$BP) currchron = 1950-1000*RSLresults$chrons[sample(1:10000,ss),i]
      if(RSLresults$BP) currchron = RSLresults$chrons[sample(1:10000,ss),i]
      currRSL = rnorm(ss,mean=RSLresults$RSLmean[i],sd=RSLresults$RSLsd[i])
      out = kde2d(currchron,currRSL,n=100)
      image(out,add=TRUE,col=COLOURS,xlim=par('usr')[1:2],ylim=par('usr')[3:4])
    }
  }
  if(ellipses) {
    for(i in 1:length(RSLresults$RSLmean)) {
      agescale = (age.high[i]-age.low[i])/4
      rslscale = RSLresults$RSLsd[i]
      lines(ellipse(0,scale=c(agescale,rslscale),centre=c(age.med[i],RSLresults$RSLmean[i])),col=rgb(0,0,1,0.4))
    }
  }
  
  if(is.null(xgrid)) {
    xgrid = seq(max(age.low),min(age.high),length=100)
  }
  
  pred.lines = matrix(NA,ncol=length(xgrid),nrow=nrow(RSLresults$samples))
  degmat = matrix(rep(0:(RSLresults$degree),length(xgrid)*(RSLresults$degree+1)),nrow=length(xgrid),ncol=RSLresults$degree+1,byrow=TRUE)
  X.pred = matrix(rep(xgrid-RSLresults$const,RSLresults$degree+1),ncol=RSLresults$degree+1)
  X.pred = X.pred^degmat
  
  for(i in 1:nrow(pred.lines)) {
    pred.lines[i,] = X.pred%*%matrix(RSLresults$samples[i,],ncol=1,nrow=RSLresults$degree+1)
  }

  pred.med = apply(pred.lines,2,'quantile',probs=0.5)
  pred.low = apply(pred.lines,2,'quantile',probs=0.025)
  pred.high = apply(pred.lines,2,'quantile',probs=0.975)
  lines(xgrid,pred.med,col='red',lwd=2)
  lines(xgrid,pred.low,col='red',lwd=2,lty=2)
  lines(xgrid,pred.high,col='red',lwd=2,lty=2)



}
