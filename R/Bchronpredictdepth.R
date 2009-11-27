Bchronpredictdepth <- function(Bchrondata,age,numchron=10000,nbreaks=30,histcolour="light blue",...) {

cat("Predict depths for a set age in the core. \n")

ndet <- nrow(Bchrondata$input)

Bchrondata$chrons <- as.matrix(read.table(Bchrondata$chronsfile))

xdepths <- rep(0,numchron)

# For each chronology, loop through and get a depth for the specified age
if(numchron>nrow(Bchrondata$chrons)) stop("Not enough chronologies; reduce numchron")
for(i in 1:numchron) {
    xdepths[i] <- approx(Bchrondata$chrons[i,],Bchrondata$outdepths,xout=age)$y
}

# create some plots
newgraphwindow(...)
hist(xdepths,main=paste(Bchrondata$fullname,": ",age," k cal years BP",sep=""),freq=FALSE,col=histcolour,xlab="Depth (cm)",las=1,breaks=nbreaks)
grid()
mtext(paste("Bchron",ifelse(Bchrondata$version>0,paste(" v",Bchrondata$version),""),sep=""),side=1,line=4,adj=1,cex=0.6)
mtext(paste(Bchrondata$calname),side=1,line=4,adj=0,cex=0.6)

}
