Bchronchangecalcurve <- function(Bchrondata) {

cat("Change user calibration curve. \n\n")
if(Bchrondata$SHOULDRUN == FALSE) {
    cat("No existing calibration curve found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

cat("Files with .bch as the file extension have already been converted for Bchron use. \n")
cat(".14c files can be downloaded directly from eg \n")
cat("http://www.radiocarbon.org/IntCal04.htm \n")
cat("but will need to be converted. This will take a few moments.\n")
choices <- c(list.files(paste(Bchrondata$path,"/CalCurve",sep="")))
title <- "List of available calibration curves:"
choose <- menu(choices, title = title)

Bchrondata$calibcurvefile <- paste(Bchrondata$path,"/CalCurve/",choices[choose],sep="")

if(substring(Bchrondata$calibcurvefile,nchar(Bchrondata$calibcurvefile)-2,nchar(Bchrondata$calibcurvefile)) == "bch")
{
cat("Calibration curve set. \n")
return(Bchrondata)
} else if(substring(Bchrondata$calibcurvefile,nchar(Bchrondata$calibcurvefile)-2,nchar(Bchrondata$calibcurvefile)) == "14c")
{

cat("Converting calibration curve... \n \n")

# Interpolate and create a big file for the two calibration curve look-ups
# 1 - Given cal years BP return 14C years BP
# 2 - Given cal years BP return error in years BP

# First read in file
CalCurveTemp <- read.table(Bchrondata$calibcurvefile,header=FALSE,sep=",")[,1:3]

BigCalBP <- seq(min(CalCurveTemp[,1]),max(CalCurveTemp[,1]),by=1)

Bchrondata$lowcal <- min(CalCurveTemp[,1])/1000
Bchrondata$highcal <- max(CalCurveTemp[,1])/1000

BigC14BP <- approx(CalCurveTemp[,1],CalCurveTemp[,2],xout=BigCalBP,rule=2)$y
BigSigBP <- approx(CalCurveTemp[,1],CalCurveTemp[,3],xout=BigCalBP,rule=2)$y

Bchrondata$calibcurvefile <- substring(Bchrondata$calibcurvefile,1,nchar(Bchrondata$calibcurvefile)-3)
Bchrondata$calibcurvefile <- paste(Bchrondata$calibcurvefile,"bch",sep="")

for(i in 1:length(BigCalBP)) {
    if(i%%100==0) print(i)
    cat(BigC14BP[i]/1000,BigSigBP[i]/1000,"\n",file=Bchrondata$calibcurvefile,append=ifelse(i==1,FALSE,TRUE))
}

cat("Calibration curve set. \n")
BigCalTemp <- read.table(Bchrondata$calibcurvefile)
Bchrondata$bigcalsize <- length(BigCalTemp[,1])
return(Bchrondata)

} else {

cat("Calibration curve not recognised, \n")
return(Bchrondata)

}

}