`PlotDens` <-
function(datesfile,depthfile,fullname,outfile,vers,ngroup)
{

# Needs the hdrcde library
library(hdrcde)

# Read in dates
dates <- read.table(paste(datesfile))

# Read in depths
if(ngroup==0) {
    depths <- read.table(paste(depthfile))
} else {
    depths <- matrix(1,nrow=1,ncol=1) 
} 

# create some plots
for(i in 1:nrow(depths))
{
    if(i!=1) {
        if(R.Version()$arch=="i386") {
            windows()
        } else {
            quartz()
        }
    }
    if(ngroup==0) hdr2.den(dates[,i],h=bw.nrd0(dates[,i]),main=paste(fullname,": depth = ",depths[i,1],"cm",sep=""),xlab="k cal yrs BP")
    if(ngroup>0) hdr2.den(dates[,i],h=bw.nrd0(dates[,i]),main=paste(fullname,": event = ",ngroup,sep=""),xlab="k cal yrs BP")
    mtext(paste("Bchron v",vers),side=1,line=4,adj=1,cex=0.6)
    
    # Depending on the skew, put the legend in the top left or top right of the plot    
    n = length(dates[,i])
    cubes = (dates[,i]-mean(dates[,i]))^3
    squares = (dates[,i]-mean(dates[,i]))^2
    
    skew = (sqrt(n)*sum(cubes))/(sum(squares)^(3/2))
    skew = ((sqrt(n*(n-1))/(n-2))/(n-2))*skew
    
    # Do a legend
    legend(ifelse(skew<0,min(dates[,i]),quantile(dates[,i],0.95)),max(density(dates[,i])$y),legend=c("99% HDR","95% HDR","50% HDR"),col=c("blue","red","green"),lwd=c(5,5,5),bty="n",text.width=0.8)

    # Output the 95% hdrs to a suitable file
    hdr = hdr(dates[,i],h=bw.nrd0(dates[,i]))$hdr[2,]
    hdr = hdr[!is.na(hdr)]
    cat(depths[i,1],hdr,"\n",file=outfile,append=ifelse(i==1,FALSE,TRUE))
}


}

