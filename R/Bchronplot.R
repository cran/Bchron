Bchronplot <- function(Bchrondata,plot=TRUE,dates=TRUE,datelabels=FALSE,chrons=TRUE,dateselect=NULL,limits=NULL,dateheight=mean(diff(Bchrondata$input[,4]))/200,chronwidth=0.95,datecolour="black",chroncolour="blue",outcolour="red",transp=0.5,outprob=0.1,legendloc="topleft",leghoriz=TRUE,...) {

if(dates==FALSE & chrons==FALSE) stop("One of dates of chrons must be set as true ")

# Get the chronology output
cat("Reading in data...\n")

if(dates==TRUE) {
    if(!file.exists(Bchrondata$calibdatesfile)) stop("Calibrated dates not found")
    Bchrondata$calibdates <- as.matrix(read.table(Bchrondata$calibdatesfile)) 
}

if(chrons==TRUE) {
    if(!file.exists(Bchrondata$chronsfile)) stop("Chronologies file not found")
    Bchrondata$chrons <- as.matrix(read.table(paste(Bchrondata$chronsfile)))
}

tmp <- match(dateselect,seq(1,nrow(Bchrondata$input)))
if(any(is.na(tmp))) stop("Selected dates do not match dates in input file")

# Get info from input file
depth <- Bchrondata$input[,4]/100
types <- as.integer(Bchrondata$input[,8])
ndet <- nrow(Bchrondata$input)
if(is.null(dateselect)) dateselect <- seq(1,ndet)

# Get desired dimensions for plotting.
# limits should be a 4-vector with dimensions as c(lowestdepth,oldestage,highestdepth,youngestage)
# in thousands of years, and m where appropriate
if(is.null(limits)) {
    if(dates==TRUE) {
        if(var(Bchrondata$calibdates[,ncol(Bchrondata$calibdates)])>0) {
            myhdrhigh <- hdr(Bchrondata$calibdates[,ncol(Bchrondata$calibdates)],h=bw.nrd0(Bchrondata$calibdates[,ncol(Bchrondata$calibdates)]),prob=chronwidth*100)$hdr[1,]
        } else {
            myhdrhigh <- Bchrondata$calibdates[1,ncol(Bchrondata$calibdates)]
        }
        if (var(Bchrondata$calibdates[,1])>0) {
            myhdrlow <- hdr(Bchrondata$calibdates[,1],h=bw.nrd0(Bchrondata$calibdates[,1]),prob=chronwidth*100)$hdr[1,]        
        } else {
            myhdrlow <- Bchrondata$calibdates[1,1]
        }
    }
    if(chrons==TRUE) {
        if(var(Bchrondata$chrons[,ncol(Bchrondata$chrons)])>0) {
            myhdrhigh <- hdr(Bchrondata$chrons[,ncol(Bchrondata$chrons)],h=bw.nrd0(Bchrondata$chrons[,ncol(Bchrondata$chrons)]),prob=chronwidth*100)$hdr[1,]
        } else {
            myhdrhigh <- Bchrondata$chrons[1,ncol(Bchrondata$chrons)]
        }
        if(var(Bchrondata$chrons[,1])>0) {
            myhdrlow <- hdr(Bchrondata$chrons[,1],h=bw.nrd0(Bchrondata$chrons[,1]),prob=chronwidth*100)$hdr[1,]
        } else {
            myhdrlow <- Bchrondata$chrons[1,1]
        }
    }
    oldageest <- max(myhdrhigh[!is.na(myhdrhigh)])
    youngageest <- min(myhdrlow[!is.na(myhdrlow)])
    limits <- c(max(Bchrondata$outdepth/100),oldageest,min(Bchrondata$outdepth/100),youngageest)
}
YLIMIT = c(limits[1],limits[3])
XLIMIT = c(limits[2],limits[4])

# Plot the Bchron output
if(plot==TRUE) {
    newgraphwindow(...)
    plot(1,1,xlim=XLIMIT,ylim=YLIMIT,xlab="k cal yrs BP",ylab="Depth (m)",type="n",main=paste(Bchrondata$fullname),bty="n",cex.axis=0.6,las=1,xaxt="n",yaxt="n")
    if(datelabels==TRUE) {
        bigstringw <- min(strwidth((Bchrondata$input[,1])))
        XLIMIT[1] <- XLIMIT[1]-bigstringw
    }
    if(!is.null(legendloc)) {
        if(legendloc=="topleft" & chrons==TRUE) {
            bigstringh <- strheight("legend")
            YLIMIT[2] <- YLIMIT[2]+bigstringh*2
        }
    }
    plot(1,1,xlim=XLIMIT,ylim=YLIMIT,xlab="k cal yrs BP",ylab="Depth (m)",type="n",main=paste(Bchrondata$fullname),bty="n",cex.axis=0.6,las=1,xaxt="n",yaxt="n")
    if(chrons==TRUE) {
         axis(1,at=pretty(Bchrondata$chrons))
    } else { 
        axis(1,at=pretty(Bchrondata$calibdates))
    }
    axis(2,at=pretty(Bchrondata$outdepths/100),las=1)
    mtext(paste("Bchron",ifelse(Bchrondata$version>0,paste(" v",Bchrondata$version),""),sep=""),side=1,line=4,adj=1,cex=0.6)
    mtext(paste(Bchrondata$calname),side=1,line=4,adj=0,cex=0.6)

    # Draw a grid and add lines at each depth
    grid()
}

########################################################

#Plot the dates

if(chrons==TRUE) {
    # Read in outlier file
    Bchrondata$outlier <- read.table(Bchrondata$outlierfile,header=TRUE)
    maxoutprobs <- apply(Bchrondata$outlier[,2:3],1,"max")
} else {
    maxoutprobs <- rep(0,ndet)
}

if(dates==TRUE)
{
    cat("Adding calibrated dates...\n")
    for(k in dateselect) {
        cat(k,"\n")
        
        datedens <- density(Bchrondata$calibdates[,k])
        if(max(datedens$x)-min(datedens$x)>0.01) polygon(datedens$x,depth[k]-datedens$y*dateheight/max(datedens$y),col=ifelse(maxoutprobs[k]>outprob,outcolour,datecolour),border=ifelse(maxoutprobs[k]>outprob,outcolour,datecolour))
        if(max(datedens$x)-min(datedens$x)<0.01) points(mean(datedens$x),depth[k],pch=20,col=datecolour)
    }
    
    for(k in dateselect) {
        if(datelabels==TRUE) {
            text(par("usr")[1],depth[k],Bchrondata$input[k,1],cex=0.8,pos=4,col=ifelse(maxoutprobs[k]>outprob,outcolour,datecolour))
            lines(c(par("usr")[1]+strwidth(Bchrondata$input[k,1]),par("usr")[2]),c(depth[k],depth[k]),col=ifelse(maxoutprobs[k]>outprob,outcolour,datecolour)) 
        } else {
            abline(h=depth[k],col=ifelse(maxoutprobs[k]>outprob,outcolour,datecolour))
        }
    }     
        
}


########################################################

# Add the chronologies
if(chrons==TRUE) {
    cat("Adding chronologies...\n")
    HPD <- matrix(0,nrow=length(Bchrondata$outdepths),ncol=2)
    tmp <- col2rgb(chroncolour)
    mycol <- rgb(tmp[1,1]/255,tmp[2,1]/255,tmp[3,1]/255)
    mycol2 <- paste(mycol,as.character(as.hexmode(round(transp*255,0))),sep="")
    for(j in 1:length(Bchrondata$outdepths)) {
        cat(j,"\n")
        if(var(Bchrondata$chrons[,j])>0) {
            currenthdr <- hdr(Bchrondata$chrons[,j],h=bw.nrd0(Bchrondata$chrons[,j]),prob=chronwidth*100)$hdr
        } else {
            currenthdr <- c(Bchrondata$chrons[,j],Bchrondata$chrons[,j])
        }
    
        HPD[j,] <- c(min(currenthdr),max(currenthdr))    
    }
    
    polygon(c(HPD[,1],rev(HPD[,2])),c(Bchrondata$outdepths/100,rev(Bchrondata$outdepths/100)),col=mycol2,border=mycol2)

    baddates <- sum(maxoutprobs>outprob)

# Add a legend
if(!is.null(legendloc)) {
    if(baddates==0) {
        legend(legendloc,c("Calibrated dates",paste(round(chronwidth*100,1),"% chronology",sep="")),pch=15,col=c(datecolour,mycol2),cex=0.8,bty="n",horiz=leghoriz)
    } else {
        legend(legendloc,c("Calibrated dates","Outlying dates",paste(round(chronwidth*100,1),"% chronology",sep="")),cex=0.8,pch=15,col=c(datecolour,outcolour,mycol2),bty="n",horiz=leghoriz)
    }
}

}

cat("Completed!\n")
}
