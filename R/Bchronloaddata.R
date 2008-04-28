Bchronloaddata <- function(Bchronversion=0) {

Bchrondata <- list()
Bchrondata$version <- Bchronversion
Bchrondata$EXIT <- FALSE
Bchrondata$SHOULDRUN <- FALSE
Bchrondata$CALIBRATED <- FALSE

cat("Load in a Bchron data file. \n")
cat("You need to have created an appropriate file as outlined in the help section. \n")
cat("\n")

# Get the path
BADPATH <- TRUE
while(BADPATH==TRUE) {
  cat("Please enter the path to the directory that contains the data files \n")
  cat("and the calibration curve, eg C:/cores. \n")
  cat("Leave blank for the default of C:/Bchron \n")
  Bchrondata$path <- scan(what = "", nlines = 1, quiet = TRUE)
  if(length(Bchrondata$path) == 0) Bchrondata$path <- "C:/Bchron"
  if(Bchrondata$path==0) return(Bchrondata)
  if(!file.exists(Bchrondata$path)) {
    cat("Error: path cannot be found. \n")
    cat(paste("Currently it is set at: ",Bchrondata$path,"\n", sep = ""))
    cat("Press <Enter> to try again...")
    readline()
    invisible()
  } else {
    BADPATH <- FALSE
  }
}

# Get calibration curve
BADCALIB <- TRUE
while(BADCALIB == TRUE) {
  cat("Please enter the name of the calibration curve .bch file you wish to use. \n")
  cat("Leave blank for the default of IntCal04. \n")
  cat("You can change the calibration curve later from the main menu. \n")
  calibname <- scan(what = "", nlines = 1, quiet = TRUE)
  if(length(calibname)==0) calibname <- "IntCal04"
  if(calibname==0) return(Bchrondata)
  Bchrondata$calibcurvefile <- paste(Bchrondata$path,"/CalCurve/",calibname,".bch",sep = "")
  if(!file.exists(Bchrondata$calibcurvefile)) {
    cat("Error: Calibration curve cannot be found. \n")
    cat("Check the path you have specified. \n")
    cat(paste("Currently it is set at: ",Bchrondata$calibcurvefile,"\n", sep = ""))
    cat("Press <Enter> to try again...")
    readline()
    invisible()
  } else {
    BADCALIB <- FALSE
    BigCalTemp <- read.table(Bchrondata$calibcurvefile)
    Bchrondata$bigcalsize <- length(BigCalTemp[,1])
    Bchrondata$lowcal <- -5/1000
    Bchrondata$highcal <- 26
    if(calibname!="IntCal04") {
      BigCalTemp2 <- read.table(paste(Bchrondata$path,"/CalCurve/",calibname,".14c",sep = ""),header=FALSE,sep=",")
      Bchrondata$highcal <- max(BigCalTemp2[,1])
      Bchrondata$lowcal <- min(BigCalTemp2[,1])
    }
    
  }
}

# Input file
BADINPUT <- TRUE
while(BADINPUT == TRUE) {
  cat("Now please enter the name of the .dat file you wish to run Bchron from. \n")
  BADNAME <- TRUE
  while(BADNAME==TRUE) {
    Bchrondata$name <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(Bchrondata$name)==0) Bchrondata$name <- scan(what = "", nlines = 1, quiet = TRUE)
    if(length(Bchrondata$name)>0) BADNAME <- FALSE
    if(Bchrondata$name==0) return(Bchrondata)
  }
  
  Bchrondata$ddepthfile <- paste(Bchrondata$path, "/Input/",Bchrondata$name,"ddepths.txt",sep="")
  Bchrondata$inputfile <- paste(Bchrondata$path, "/Input/",Bchrondata$name,".dat",sep="")
  if(!file.exists(Bchrondata$inputfile) || !file.exists(Bchrondata$ddepthfile)) {
      cat("Error: cannot find input/design depth files. \n")
      cat(paste("Check both ",Bchrondata$inputfile," and ",Bchrondata$ddepthfile,"\n", sep = ""))
      cat("exist and can be read.\n")
      cat("Press <Enter> to try again...")
      readline()
      invisible()
  } else {
    BADINPUT <- FALSE
  }
}

# Fullname
cat("Please enter the full name of the core (used for graph titles; leave blank if the same as above). \n")
Bchrondata$fullname <- scan(what = "", nlines = 1, quiet = TRUE,sep = "\t")
if(length(Bchrondata$fullname) == 0) Bchrondata$fullname <- Bchrondata$name
if(Bchrondata$fullname==0) return(Bchrondata)
cat("Thank you.\n\n")

cat("Data entered correctly. \n")
cat("Press <Enter> to continue...")
readline()
invisible()

Bchrondata$calibdatesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"TrueDates.txt",sep = "")
Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "pars.txt",sep = "")
if(file.exists(Bchrondata$calibdatesfile)) Bchrondata$CALIBRATED <- TRUE
if(file.exists(Bchrondata$parsfile)) Bchrondata$RUN <- TRUE
Bchrondata$SHOULDRUN <- TRUE
return(Bchrondata)

}