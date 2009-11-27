Bchroncreatecalcurve <- function() {

cat("Create user calibration curve. \n\n")

BADPATH <- TRUE
while(BADPATH==TRUE) {
  cat("First, please enter the path to the Bchron directory that contains the data files \n")
  cat("and the calibration curve. \n")
  if(.Platform$OS.type=="unix") {
      cat("Leave blank for the default of",paste(getwd(),"/Bchron",sep=""),"\n")
  } else {
      cat("Leave blank for the default of C:/Bchron \n")
  }
  path <- scan(what = "", nlines = 1, quiet = TRUE)
  if(length(path) == 0) {
    if(.Platform$OS.type=="unix") {
        path <- paste(getwd(),"/Bchron",sep="")
    } else {
        path <- "C:/Bchron"
    }
  }
  if(!file.exists(path)) {
    cat("Error: path cannot be found. \n")
    cat(paste("Currently it is set at: ",path,"\n", sep = ""))
    cat("Press <Enter> to try again...")
    readline()
    invisible()
  } else {
    BADPATH <- FALSE
  }
}

cat("Files with .bch as the file extension have already been converted for Bchron use. \n")
cat(".14c files can be downloaded directly from eg \n")
cat("http://www.radiocarbon.org \n")
cat("but will need to be converted. This will take a few moments.\n")
choices <- c(list.files(paste(path,"/CalCurve",sep="")))
title <- "List of available calibration curves:"
choose <- menu(choices, title = title)

calibcurvefile <- paste(path,"/CalCurve/",choices[choose],sep="")

if(substring(calibcurvefile,nchar(calibcurvefile)-2,nchar(calibcurvefile)) == "bch") {
    cat("Calibration curve already created. \n")
} 
else if(substring(calibcurvefile,nchar(calibcurvefile)-2,nchar(calibcurvefile)) == "14c") {

cat("Converting calibration curve... \n \n")

# Interpolate and create a big file for the two calibration curve look-ups
# 1 - Given cal years BP return 14C years BP
# 2 - Given cal years BP return error in years BP

# First read in file
CalCurveTemp <- read.table(calibcurvefile,header=FALSE,sep=",")[,1:3]

BigCalBP <- seq(min(CalCurveTemp[,1]),max(CalCurveTemp[,1]),by=1)
BigC14BP <- approx(CalCurveTemp[,1],CalCurveTemp[,2],xout=BigCalBP,rule=2)$y
BigSigBP <- approx(CalCurveTemp[,1],CalCurveTemp[,3],xout=BigCalBP,rule=2)$y

calibcurvefile <- substring(calibcurvefile,1,nchar(calibcurvefile)-3)
calibcurvefile <- paste(calibcurvefile,"bch",sep="")

for(i in 1:length(BigCalBP)) {
    if(i%%100==0) print(i)
    cat(BigC14BP[i]/1000,BigSigBP[i]/1000,"\n",file=calibcurvefile,append=ifelse(i==1,FALSE,TRUE))
}

# Plot it in case there are any inconsistencies
plot(CalCurveTemp[,1],CalCurveTemp[,2],xlab="Cal BP",ylab="14C BP",main=calibcurvefile,type="l",col="blue",las=1)
grid()
lines(CalCurveTemp[,1],CalCurveTemp[,2]-2*CalCurveTemp[,3],col="blue",lty=2)
lines(CalCurveTemp[,1],CalCurveTemp[,2]+2*CalCurveTemp[,3],col="blue",lty=2)

cat("Calibration curve set. \n")

} else {

    cat("Calibration curve not recognised, \n")

}

}
