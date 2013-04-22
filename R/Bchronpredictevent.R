Bchronpredictevent <-
function(Bchrondata,event="Event",numchron=10000) {

cat("Predict ages for an event in the core. \n")

# Get the number of determinations
ndet <- nrow(Bchrondata$input)

eventdepthfile <- paste(Bchrondata$path,"/Input/",Bchrondata$name,"EventDepths",event,".txt",sep="")
if(!file.exists(eventdepthfile)) stop(paste("Event depths file not found: ",eventdepthfile,"\n"))
eventagefile <- paste(Bchrondata$path,"/Output/",Bchrondata$name,"EventAges",event,".txt",sep="")
evdepth <- read.table(eventdepthfile)/100

lowdepths <- evdepth[, 1]
highdepths <- evdepth[, 2]
nddepthints <- length(lowdepths)
  
cat("Parameters file at", Bchrondata$parsfile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Event depths file at", eventdepthfile, " \n")
cat("Event ages file created at", eventagefile, " \n")
cat("Range of output depths is \n")
for(j in 1:nddepthints) {
    cat(lowdepths[j],"to",highdepths[j],"\n")
}
cat("Outlier output file at", Bchrondata$outlierfile, " \n")
cat("Number of determinations is", ndet, " \n")
cat("Number of event depth intervals is", nddepthints, " \n")
cat("Number of samples is", numchron, " \n")
cat("Date of extraction of core (in k yrs BP) is",Bchrondata$extractdate, " \n")

cat("\n")

out <- .C("predictrand",
    as.character(Bchrondata$input[,1]),
    as.double(Bchrondata$input[,2]/1000),
    as.double(Bchrondata$input[,3]/1000),
    as.double(Bchrondata$input[,4]/100),
    as.double(Bchrondata$input[,5]/100),
    as.double(Bchrondata$input[,6]),
    as.double(Bchrondata$input[,7]),
    as.integer(Bchrondata$input[,8]),
    as.character(Bchrondata$parsfile),
    as.character(eventagefile),
    as.double(lowdepths),
    as.double(highdepths),
    as.integer(nddepthints),
    as.integer(ndet),
    as.integer(numchron),
    as.double(Bchrondata$extractdate),
    as.character(Bchrondata$outlierfile))

cat("Completed!\n")

}
