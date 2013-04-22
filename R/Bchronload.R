Bchronload <-
function(name,fullname=NULL,path=NULL,outdepths=NULL,calibname="IntCal09",extractdate=-0.05,check=FALSE,full=FALSE) {

cat("Loading data...\n")

# Change default path if on a mac
if(is.null(path)) {
   if(.Platform$OS.type=="unix") { path <- paste(getwd(),"/Bchron",sep="") }
    else { path <- "C:/Bchron" }
}

# Create an object and fill in the components
Bchrondata <- list()

# Get the various names
Bchrondata$name <- name
Bchrondata$fullname <- fullname
if(is.null(fullname)) Bchrondata$fullname <- Bchrondata$name
Bchrondata$version <- try(read.dcf(file=system.file("DESCRIPTION", package="Bchron"),fields="Version"),silent=TRUE)

# Find the various files
Bchrondata$path <- path
Bchrondata$calname <- calibname
Bchrondata$calibcurvefile <- paste(Bchrondata$path,"/CalCurve/",calibname,".bch",sep = "")
Bchrondata$c14file <- paste(Bchrondata$path,"/CalCurve/",calibname,".14c",sep = "")
Bchrondata$inputfile <- paste(Bchrondata$path, "/Input/",Bchrondata$name,".dat",sep="")
Bchrondata$calibdatesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"TrueDates.txt",sep = "")
Bchrondata$calibrangesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"CalibratedRanges.txt",sep = "")
Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "pars.txt",sep = "")
Bchrondata$chronsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "chrons.txt",sep = "")
Bchrondata$outlierfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"outliers.txt", sep = "")
Bchrondata$rangesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"ranges.txt", sep = "")
Bchrondata$extractdate <- extractdate

# Get some details regarding the calibration curve dimensions
BigCalTemp <- read.table(Bchrondata$calibcurvefile)
BigC14Temp <- read.table(Bchrondata$c14file,sep=",")
Bchrondata$bigcalsize <- length(BigCalTemp[,1])
Bchrondata$lowcal <- min(BigC14Temp[,1])/1000
Bchrondata$highcal <- max(BigC14Temp[,1])/1000

Bchrondata$input <- read.table(Bchrondata$inputfile,header=TRUE)
Bchrondata$input[,1] <- as.character(Bchrondata$input[,1])
if(is.null(outdepths)) {
    cat("Using default range of top depth to bottom depth for output depths. \n")
    outdepths <- seq(min(Bchrondata$input[,4]),max(Bchrondata$input[,4]),length=200)
}
Bchrondata$outdepths <- outdepths
if(length(outdepths)>400) cat("WARNING: greater than 400 output depths given; may result in slow prediction runs. \n")

if(full==TRUE) {
    cat("Loading in parameters, calibrated dates, chronologies, outliers, and ranges (where available) \n")
    if(file.exists(Bchrondata$calibdatesfile)) Bchrondata$calibdates <- as.matrix(read.table(Bchrondata$calibdatesfile))
    if(file.exists(Bchrondata$calibrangesfile)) Bchrondata$calibranges <- read.table(Bchrondata$calibrangesfile,header=TRUE)
    if(file.exists(Bchrondata$parsfile)) Bchrondata$pars <- as.matrix(read.table(Bchrondata$parsfile))
    if(file.exists(Bchrondata$chronsfile)) Bchrondata$chrons <- as.matrix(read.table(Bchrondata$chronsfile))
    if(file.exists(Bchrondata$outlierfile)) Bchrondata$outlier <- read.table(Bchrondata$outlierfile,header=TRUE)
    if(file.exists(Bchrondata$rangesfile)) Bchrondata$ranges <- read.table(Bchrondata$rangesfile,header=TRUE)
}

cat("Data loaded. \n")
cat("Note: path to Bchron set as:",path,"\n")
cat("\n")

if(check==TRUE) {
    cat("Now checking files for possible errors. \n")
    cat("Press <Enter> to continue or <Esc> to exit...")
    readline()
    invisible()
    Bchroncheck(Bchrondata)
}

class(Bchrondata) = "Bchron"
return(Bchrondata)

}
