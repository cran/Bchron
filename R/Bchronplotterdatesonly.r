Bchronplotterdatesonly <- function(Bchrondata,asklimits=FALSE) {

library(hdrcde)

cat("Plotting calibrated dates... \n")
if (Bchrondata$SHOULDRUN == FALSE || Bchrondata$CALIBRATED == FALSE) {
    cat("No data found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

Bchrondata$datesfile <- paste(Bchrondata$path,"/Output/",Bchrondata$name, "TrueDates.txt",sep = "")

choices3 <- c("colour", "black and white")
choose3 <- menu(choices3, title = "How would you like the dates to be drawn?")
Bchrondata$COLOUR <- ifelse(choose3 == 1, TRUE, FALSE)

# Get the chronology output
cat("Reading in data...\n")

# Get the unconstrained dates
Bchrondata$calibdatesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"TrueDates.txt",sep = "")
TrueDates <- as.matrix(read.table(Bchrondata$calibdatesfile))

# Get the design depths
ddepths <- read.table(Bchrondata$ddepthfile,header=FALSE)/100

# Get the depths
temp <- read.table(Bchrondata$inputfile,header=TRUE)
depth <- temp[,4]/100

# Now do ages
myhdr <- hdr(TrueDates[,ncol(TrueDates)],h=bw.nrd0(TrueDates[,ncol(TrueDates)]))$hdr[1,]
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

#Plot TrueDates
# Draw some graphs of thetas too.
BigThetaAll <- TrueDates
cat("Adding calibrated dates...\n")
for(k in 1:length(depth)) {
    print(k)

    #if(depth[k]!=0){
      HDR <- try(hdr(BigThetaAll[,k],h=bw.nrd0(BigThetaAll[,k]))$hdr[2,],silent=TRUE)

      if(length(HDR)>1) {
        HDR <- HDR[!is.na(HDR)]
        HDRcol <- length(HDR)

        if(Bchrondata$COLOUR==TRUE) for(t in seq(1,HDRcol,by=2)) lines(c(HDR[t],HDR[t+1]),c(depth[k],depth[k]),col="blue",lwd=6)
        if(Bchrondata$COLOUR==FALSE) for(t in seq(1,HDRcol,by=2)) lines(c(HDR[t],HDR[t+1]),c(depth[k],depth[k]),col="grey15",lwd=6)
      }
    #}
}
#if(Bchrondata$COLOUR==TRUE) legend(x=XLIMIT[1],y=YLIMIT[2],legend=c("95% HDR","75% HDR","50% HDR","Mode","Unrestricted Dates"),col=c(colours[1],colours[2],colours[3],"black","blue"),lwd=c(lwidths[1],lwidths[2],lwidths[3],-1,6),lty=c(1,1,1,-1,1),pch=c(-1,-1,-1,19,-1),bty="n",text.width=0.8)
#if(Bchrondata$COLOUR==FALSE) legend(x=XLIMIT[1],y=YLIMIT[2],legend=c("95% HDR","Mode","Unrestricted Dates"),col=c("lightgrey","black","grey15"),lwd=c(3,-1,6),lty=c(1,-1,1),pch=c(-1,19,-1),bty="n",text.width=0.8)

  
cat("Completed!\n")
cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

return(Bchrondata)
}
