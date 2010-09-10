Bchronproxyplot <-
function(Bchrondata,proxy,title=NULL,xlabel="Age (k cal yrs BP)",ylabel="Proxy",num=10,col=c("red","red","red"),smooth=FALSE,...){

# Read in all chronologies  
Bchrondata$chrons <- as.matrix(read.table(Bchrondata$chronsfile))
Bchrondata$ranges <- read.table(Bchrondata$rangesfile,header=TRUE)

# Now produce plot of lines
xlimits <- c(min(Bchrondata$ranges[,2]),max(Bchrondata$ranges[,4]))
ylimits <- c(min(proxy),max(proxy))
dev.new(...)
plot(1,1,type="n",xlim=xlimits,ylim=ylimits,main=title,xlab= xlabel,ylab= ylabel,las=1)
grid()

# Plot num lines as well
myrows <- sample(seq(1,nrow(Bchrondata$chrons)),num,replace=FALSE)
if(num>0) for(i in 1:num) lines(Bchrondata$chrons[myrows[i],],proxy,col="grey")

# Plot an estimate of the median line
if(smooth==FALSE) {
lines(Bchrondata$ranges[,2],proxy,col=col[1],lty="dotted",lwd=2)
lines(Bchrondata$ranges[,3],proxy,col=col[2],lwd=2)
lines(Bchrondata$ranges[,4],proxy,col=col[3],lty="dotted",lwd=2)
} else {
tmp <- ksmooth(Bchrondata$ranges[,2],proxy,bandwidth=bw.nrd0(Bchrondata$ranges[,2]))
tmp2 <- ksmooth(Bchrondata$ranges[,3],proxy,bandwidth=bw.nrd0(Bchrondata$ranges[,3]))
tmp3 <- ksmooth(Bchrondata$ranges[,4],proxy,bandwidth=bw.nrd0(Bchrondata$ranges[,4]))
lines(tmp,col=col[1],lty="dotted",lwd=2)
lines(tmp2,col=col[2],lwd=2)
lines(tmp3,col=col[3],lty="dotted",lwd=2)
}

}

