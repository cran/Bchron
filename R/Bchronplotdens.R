Bchronplotdens <- function(datesfile,fullname,outfile,eventname,vers=0)
{

dates <- read.table(datesfile,header=FALSE)
# create some plots
newgraphwindow()
hdr2.den(dates[,1],h=bw.nrd0(dates[,1]),main=paste(fullname,": ",eventname,sep=""),xlab="k cal yrs BP")
mtext(paste("Bchron",ifelse(vers>0,paste(" v",vers),""),sep=""),side=1,line=4,adj=1,cex=0.6)

# Depending on the skew, put the legend in the top left or top right of the plot    
n = length(dates[,1])
cubes = (dates-mean(dates[,1]))^3
squares = (dates-mean(dates[,1]))^2
skew = (sqrt(n)*sum(cubes))/(sum(squares)^(3/2))
skew = ((sqrt(n*(n-1))/(n-2))/(n-2))*skew

# Do a legend
legend(ifelse(skew<0,min(dates[,1]),quantile(dates[,1],0.95)),max(density(dates[,1])$y),legend=c("99% HDR","95% HDR","50% HDR"),col=c("blue","red","green"),lwd=c(5,5,5),bty="n",text.width=0.8)

# Output the 95% hdrs to a suitable file
hdr = hdr(dates[,1],h=bw.nrd0(dates[,1]))$hdr[2,]
hdr = hdr[!is.na(hdr)]*1000
cat(paste("95% HDR for event:",eventname,"\n"))
cat(hdr,"\n")
cat(hdr,"\n",file=outfile,append=FALSE)

}

