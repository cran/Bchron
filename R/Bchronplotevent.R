Bchronplotevent <-
function(Bchrondata,depth=NULL,slice=NULL,event=NULL,eventname=NULL,nbreaks=30,histcolour="light blue",...)
{

if(is.null(depth) & is.null(slice) & is.null(event)) {
    choices <- c("slice","depth", "event")
    choose <- menu(choices,title="Please choose which of 'slice', 'depth' or 'event' you wish to enter")
    
    if(choose==1) {
        cat("Please enter slice number. Choose an integer from 1 to",length(Bchrondata$outdepths),"\n")
        slice <- scan(what = "", nlines = 1, quiet = TRUE)
        while(length(slice)==0 | !is.integer(slice) | slice<1 | slice>length(Bchrondata$outdepths)) {
            cat("Invalid value, try again \n")
            slice <- scan(what = "", nlines = 1, quiet = TRUE)
        }
    }
    if(choose==2) {
        cat("Please enter depth in cm. Choose a number from the set of output depth values given when loading the data in \n")
        depth <- scan(what = "", nlines = 1, quiet = TRUE)
        while(length(depth)==0 | is.na(match(depth,Bchrondata$outdepths))) {
            cat("Invalid value, try again \n")
            depth <- scan(what = "", nlines = 1, quiet = TRUE)
        }
    }
    if(choose==3) {
        cat("Please enter event name. \n")
        event <- scan(what = "", nlines = 1, quiet = TRUE)
        while(length(event)==0) event <- scan(what = "", nlines = 1, quiet = TRUE)
        event <- as.character(event)
    }
}

if(!is.null(depth) & !is.null(slice) & !is.null(event)) stop("Only one of depth, slice, event can be given")

if(is.null(slice) & is.null(event)) {
    slice <- match(depth,Bchrondata$outdepths)
    if(is.na(slice)) stop(paste("Depth given not found output depths."))
}
if(is.null(depth) & is.null(event)) {
    depth <- Bchrondata$outdepths[slice]
}

# If given a depth/slice
if(!is.null(slice)) {

    # Get the chronologies
    Bchrondata$chrons <- as.matrix(read.table(Bchrondata$chronsfile))

    # create some plots
    dev.new(...)
    hist(Bchrondata$chrons[,slice],main=paste(Bchrondata$fullname,": ",depth," cm",sep=""),freq=FALSE,col=histcolour,xlab="k cal yrs BP",las=1,breaks=nbreaks)
    grid()
    mtext(paste("Bchron",ifelse(Bchrondata$version>0,paste(" v",Bchrondata$version),""),sep=""),side=1,line=4,adj=1,cex=0.6)
    mtext(paste(Bchrondata$calname),side=1,line=4,adj=0,cex=0.6)
}

if(!is.null(event)) {
    if(is.null(eventname)) eventname <- event
    
    eventagefile <- paste(Bchrondata$path,"/Output/",Bchrondata$name,"EventAges",event,".txt",sep="")
    if(!file.exists(eventagefile)) stop(paste("Event age file not found:",eventagefile))
    eventage <- as.matrix(read.table(eventagefile))

    # create some plots
    dev.new(...)
    hist(eventage,main=paste(Bchrondata$fullname,": ",eventname,sep=""),freq=FALSE,col=histcolour,xlab="k cal yrs BP",las=1,breaks=nbreaks)
    grid()
    mtext(paste("Bchron",ifelse(Bchrondata$version>0,paste(" v",Bchrondata$version),""),sep=""),side=1,line=4,adj=1,cex=0.6)
    mtext(paste(Bchrondata$calname),side=1,line=4,adj=0,cex=0.6)
}


}

