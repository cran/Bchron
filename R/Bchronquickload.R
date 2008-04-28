Bchronquickload <- function(name,fullname=NULL,path="c:/Bchron",calibname="IntCal04",Bchronversion=0,check=FALSE) {

Bchrondata <- list()
Bchrondata$version <- Bchronversion
Bchrondata$path <- path
Bchrondata$calibcurvefile <- paste(Bchrondata$path,"/CalCurve/",calibname,".bch",sep = "")
Bchrondata$name <- name
Bchrondata$inputfile <- paste(Bchrondata$path, "/Input/",Bchrondata$name,".dat",sep="")
Bchrondata$ddepthfile <- paste(Bchrondata$path, "/Input/",Bchrondata$name,"ddepths.txt",sep="")
if(is.null(fullname)) Bchrondata$fullname <- Bchrondata$name
Bchrondata$SHOULDRUN <- TRUE
Bchrondata$CALIBRATED <- FALSE
Bchrondata$RUN <- FALSE
Bchrondata$calibdatesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"TrueDates.txt",sep = "")
Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "pars.txt",sep = "")
if(file.exists(Bchrondata$calibdatesfile)) Bchrondata$CALIBRATED <- TRUE
if(file.exists(Bchrondata$parsfile)) Bchrondata$RUN <- TRUE
BigCalTemp <- read.table(Bchrondata$calibcurvefile)
Bchrondata$bigcalsize <- length(BigCalTemp[,1])
Bchrondata$lowcal <- -5/1000
Bchrondata$highcal <- 26
if(calibname!="IntCal04") {
  BigCalTemp2 <- read.table(paste(Bchrondata$path,"/CalCurve/",calibname,".14c",sep = ""),header=FALSE,sep=",")
  Bchrondata$highcal <- max(BigCalTemp2[,1])/1000
  Bchrondata$lowcal <- min(BigCalTemp2[,1])
}
if(check==TRUE) Bchrondata <- Bchroncheck(Bchrondata)

return(Bchrondata)

}