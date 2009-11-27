Bchronproxyplot <- function(Bchrondata,proxy,title=NULL,xlabel="Age (k cal yrs BP)",ylabel="Proxy",num=10,...){

# Read in all chronologies  
Bchrondata$chrons <- as.matrix(read.table(Bchrondata$chronsfile))
Bchrondata$ranges <- read.table(Bchrondata$rangesfile,header=TRUE)

# Now produce plot of lines
xlimits <- c(min(Bchrondata$ranges[,2]),max(Bchrondata$ranges[,4]))
ylimits <- c(min(proxy),max(proxy))
newgraphwindow(...)
plot(1,1,type="n",xlim=xlimits,ylim=ylimits,main=title,xlab= xlabel,ylab= ylabel,las=1)
grid()

# Plot num lines as well
myrows <- sample(seq(1,nrow(Bchrondata$chrons)),num,replace=FALSE)
for(i in 1:num) lines(Bchrondata$chrons[myrows[i],],proxy,col="grey")

# Plot an estimate of the median line
lines(Bchrondata$ranges[,2],proxy,col="red",lty="dotted",lwd=2)
lines(Bchrondata$ranges[,3],proxy,col="red",lwd=2)
lines(Bchrondata$ranges[,4],proxy,col="red",lty="dotted",lwd=2)

}
