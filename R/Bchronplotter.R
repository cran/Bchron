Bchronplotter <- function(Bchrondata,colours=c("yellow","orange","red"),lwidths=c(3,3,3),asklimits=FALSE) {

library(hdrcde)

cat("Plotting chronology... \n")
if (Bchrondata$SHOULDRUN == FALSE) {
    cat("No data found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

Bchrondata$datesfile <- paste(Bchrondata$path,"/Output/",Bchrondata$name, "TrueDates.txt",sep = "")
Bchrondata$rangesfile <- paste(Bchrondata$path,"/Output/",Bchrondata$name, "Ranges.txt",sep = "")
Bchrondata$chrons <- paste(Bchrondata$path,"/Output/",Bchrondata$name, "chrons.txt",sep = "")

choices3 <- c("colour", "black and white")
choose3 <- menu(choices3, title = "How would you like the chronology to be drawn?")
Bchrondata$COLOUR <- ifelse(choose3 == 1, TRUE, FALSE)

# Get the chronology output
cat("Reading in data...\n")
PredOut <- as.matrix(read.table(paste(Bchrondata$chrons)))

# Get the unconstrained dates if they exist
if(Bchrondata$CALIBRATED==TRUE) {
  Bchrondata$calibdatesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"TrueDates.txt",sep = "")
  TrueDates <- as.matrix(read.table(Bchrondata$calibdatesfile))
}

# Get the design depths
ddepths <- read.table(Bchrondata$ddepthfile,header=FALSE)/100

# Get the depths
temp <- read.table(Bchrondata$inputfile,header=TRUE)
depth <- temp[,4]/100

# Now do ages
myhdr <- hdr(PredOut[,ncol(PredOut)],h=bw.nrd0(PredOut[,ncol(PredOut)]))$hdr[1,]
topleft <- max(myhdr[!is.na(myhdr)])

# Get desired dimensions for plotting.
if(asklimits==TRUE) {
cat("Please enter depth ranges for plot (leave either of these blank for defaults) \n")
cat("Top depth of interest (in cm): \n")
topdepth <- scan(nlines=1,quiet=TRUE)
if(length(topdepth)==0) {
  topdepth <- min(ddepths)
} else {
  topdepth <- topdepth/100
}
cat("Bottom depth of interest (in cm): \n")
botdepth <- scan(nlines=1,quiet=TRUE)
if(length(botdepth)==0) {
  botdepth <- max(ddepths)
} else {
  botdepth <- botdepth/100
}

cat("Please enter age ranges for plot (again, leave either of these blank for defaults) \n")
cat("Youngest age of interest (in k cal yrs BP): \n")
youngage <- scan(nlines=1,quiet=TRUE)
if(length(youngage)==0) youngage <- 0
cat("Oldest age of interest (in k cal yrs BP): \n")
oldage <- scan(nlines=1,quiet=TRUE)
if(length(oldage)==0) oldage <- topleft
} else {
  topdepth <- min(ddepths)
  botdepth <- max(ddepths)
  youngage <- 0
  oldage <- topleft
}
YLIMIT = c(botdepth,topdepth)
XLIMIT = c(oldage,youngage)

# Plot the Bchron output
newgraphwindow()
plot(1,1,xlim=XLIMIT,ylim=YLIMIT,xlab="k cal yrs BP",ylab="Depth (m)",type="n",main=paste(Bchrondata$fullname),yaxp=c(0,30,30),xaxp=c(50,-5,55),cex.axis=0.6)
mtext(paste("Bchron",ifelse(Bchrondata$version>0,paste(" v",Bchrondata$version),""),sep=""),side=1,line=4,adj=1,cex=0.6)

# Draw some extra tick marks
axis(1,at=seq(-100,100,by=0.5),labels=FALSE,tcl=0.2)
axis(2,at=seq(-100,100,by=0.5),labels=FALSE,tcl=0.2)

# Draw (feint) gridlines
vertlines <- seq(0.0,topleft,by=1)
horlines <- seq(0,max(ddepths),by=0.5)
for(i in 1:length(vertlines)) lines(c(vertlines[i],vertlines[i]),c(-200,200),col="light grey")
for(i in 1:length(horlines)) lines(c(-200,200),c(horlines[i],horlines[i]),col="light grey")

size <- 50
HPD95 <- matrix(0,nrow=nrow(ddepths),ncol=size)
HPD75 <- matrix(0,nrow=nrow(ddepths),ncol=size)
HPD50 <- matrix(0,nrow=nrow(ddepths),ncol=size)
HPD2 <- matrix(0,nrow=nrow(ddepths),ncol=size)

cat("Adding ddepth layers...\n")
for(j in 1:nrow(ddepths)) {
    print(j)
    if(var(PredOut[,j])>0) {
        currenthdr <- hdr(PredOut[,j],h=bw.nrd0(PredOut[,j]),prob=c(50,75,95,2))$hdr
        Currenthdr95 <- currenthdr[1,]
        Currenthdr75 <- currenthdr[2,]
        Currenthdr50 <- currenthdr[3,]
        Currenthdr2 <- currenthdr[4,]
    } else {
        Currenthdr95 <- c(PredOut[,j],PredOut[,j])
        Currenthdr75 <- c(PredOut[,j],PredOut[,j])
        Currenthdr50 <- c(PredOut[,j],PredOut[,j])
        Currenthdr2 <- c(PredOut[,j],PredOut[,j])
    }


    # Get the appropriate number of columns (max 30) to fill in
    cols95 <- min(length(Currenthdr95),size)
    cols75 <- min(length(Currenthdr75),size)
    cols50 <- min(length(Currenthdr50),size)
    cols2 <- min(length(Currenthdr2),size)

    # Now put it in the matrix
    HPD95[j,1:cols95] <- Currenthdr95[1:cols95]
    HPD75[j,1:cols75] <- Currenthdr75[1:cols75]
    HPD50[j,1:cols50] <- Currenthdr50[1:cols50]
    HPD2[j,1:cols2] <- Currenthdr2[1:cols2]


    if(Bchrondata$COLOUR==TRUE) {
    # And draw some pretty lines
    for(k in 1:(cols95/2)) {
        lines(c(HPD95[j,2*k-1],HPD95[j,2*k]),c(ddepths[j,1],ddepths[j,1]),col=colours[1],lwd=lwidths[1])
    }
    for(k in 1:(cols75/2)) {
        lines(c(HPD75[j,2*k-1],HPD75[j,2*k]),c(ddepths[j,1],ddepths[j,1]),col=colours[2],lwd=lwidths[2])
    }
    for(k in 1:(cols50/2)) {
        lines(c(HPD50[j,2*k-1],HPD50[j,2*k]),c(ddepths[j,1],ddepths[j,1]),col=colours[3],lwd=lwidths[3])
    }
    for(k in 1:(cols2/2)) {
        lines(c(HPD2[j,2*k-1],HPD2[j,2*k]),c(ddepths[j,1],ddepths[j,1]),col="black",lwd=3)
    }
    }

    if(Bchrondata$COLOUR == FALSE) {
    for(k in 1:(cols95/2)) {
        lines(c(HPD95[j,2*k-1],HPD95[j,2*k]),c(ddepths[j,1],ddepths[j,1]),col="lightgrey",lwd=3)
    }
    for(k in 1:(cols2/2)) {
        lines(c(HPD2[j,2*k-1],HPD2[j,2*k]),c(ddepths[j,1],ddepths[j,1]),col="black",lwd=3)
    }
    }

    # Put some of this kind of stuff in a file so that we've got depth, 2.5%, 97.5% in three columns
    if(length(Bchrondata$rangesfile)>0) cat(ddepths[j,1],quantile(PredOut[,j]*1000,c(0.025,0.5,0.975)),"\n",file=paste(Bchrondata$rangesfile),append=ifelse(j==1,FALSE,TRUE))

}


#Plot TrueDates
# Draw some graphs of thetas too.
if(Bchrondata$CALIBRATED==TRUE)
{
BigThetaAll <- TrueDates
cat("Adding calibrated dates...\n")
for(k in 1:length(depth)) {
    print(k)

    if(depth[k]!=0){
    HDR <- try(hdr(BigThetaAll[,k],h=bw.nrd0(BigThetaAll[,k]))$hdr[2,],silent=TRUE)

    if(length(HDR)>1)
    {
     HDR <- HDR[!is.na(HDR)]
     HDRcol <- length(HDR)

    if(Bchrondata$COLOUR==TRUE) for(t in seq(1,HDRcol,by=2)) lines(c(HDR[t],HDR[t+1]),c(depth[k],depth[k]),col="blue",lwd=6)
    if(Bchrondata$COLOUR==FALSE) for(t in seq(1,HDRcol,by=2)) lines(c(HDR[t],HDR[t+1]),c(depth[k],depth[k]),col="grey15",lwd=6)
    }
    }
}
    if(Bchrondata$COLOUR==TRUE) legend(x=XLIMIT[1],y=YLIMIT[2],legend=c("95% HDR","75% HDR","50% HDR","Mode","Unrestricted Dates"),col=c(colours[1],colours[2],colours[3],"black","blue"),lwd=c(lwidths[1],lwidths[2],lwidths[3],-1,6),lty=c(1,1,1,-1,1),pch=c(-1,-1,-1,19,-1),bty="n",text.width=0.8)
    if(Bchrondata$COLOUR==FALSE) legend(x=XLIMIT[1],y=YLIMIT[2],legend=c("95% HDR","Mode","Unrestricted Dates"),col=c("lightgrey","black","grey15"),lwd=c(3,-1,6),lty=c(1,-1,1),pch=c(-1,19,-1),bty="n",text.width=0.8)

} else {
    if(Bchrondata$COLOUR==TRUE) legend(x=XLIMIT[1],y=YLIMIT[2],legend=c("95% HDR","75% HDR","50% HDR","Mode"),col=c(colours[1],colours[2],colours[3],"black"),lwd=c(lwidths[1],lwidths[2],lwidths[3],-1),lty=c(1,1,1,-1),pch=c(-1,-1,-1,19),bty="n",text.width=0.8)
    if(Bchrondata$COLOUR==FALSE) legend(x=XLIMIT[1],y=YLIMIT[2],legend=c("95% HDR","Mode"),col=c("lightgrey","black"),lwd=c(3,-1),lty=c(1,-1),pch=c(-1,19),bty="n",text.width=0.8)
}


  
cat("Completed!\n")
cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

return(Bchrondata)
}
