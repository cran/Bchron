Bchroncalibrate <- function(Bchrondata,iterations=500000,burnin=50000,thinby=45,howmany=50000) {

cat("Calibrating radiocarbon dates... \n")
if (Bchrondata$SHOULDRUN == FALSE) {
    cat("You have not entered any data. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

# Get the number of determinations
Temp <- read.table(Bchrondata$inputfile,header = TRUE)
Bchrondata$ndet <- nrow(Temp)

# Create output file path
Bchrondata$calibdatesfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"TrueDates.txt",sep = "")
if(file.exists(Bchrondata$calibdatesfile)) {
  cat("Calibration stage appears to have been already run for this core. \n")
  cat("Do you wish to re-run? (y/n) \n")
  rerun <- scan(what = "", nlines = 1, quiet = TRUE)
  while(length(rerun)==0) rerun <- scan(what = "", nlines = 1, quiet = TRUE)
  if(rerun=="n" || rerun=="no") return(Bchrondata)
}

cat("============================================================================\n")
cat("WARNING: running calibration here will over-write previous calibration runs. \n")
cat("============================================================================\n")
cat("Select 0 to exit this menu. \n")
choices2 <- c("standard", "long", "super-long")
choose2 <- menu(choices2, title = "What type of Bchron model run would you like?")
iterations <- 500000
burnin <- 50000
howmany <- 50000
thinby <- 45
if(choose2 == 2) {
    iterations <- 1e+06
    burnin <- 200000
    howmany <- 50000
    thinby <- 80
}
if (choose2 == 3) {
    iterations <- 2e+06
    burnin <- 400000
    howmany <- 50000
    thinby <- 160
}
if(choose2 == 0) return(Bchrondata)

cat("Calibration curve at", Bchrondata$calibcurvefile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Output file at", Bchrondata$calibdatesfile, " \n")
cat("Number of determinations is", Bchrondata$ndet, " \n")
cat("\n")
out <- .C("calibrate", as.character(Bchrondata$calibcurvefile),
        as.character(Bchrondata$inputfile),
        as.character(Bchrondata$calibdatesfile),
        as.integer(Bchrondata$ndet),
        as.integer(Bchrondata$bigcalsize),
        as.double(Bchrondata$lowcal),
        as.double(Bchrondata$highcal),
        as.integer(iterations),
        as.integer(burnin),
        as.integer(thinby),
        as.integer(howmany)
        )

cat("\n")
cat("Calibration completed successfully. \n")
cat(paste("Results stored in:",Bchrondata$calibdatesfile),sep="")
cat("\n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

Bchrondata$CALIBRATED <- TRUE
return(Bchrondata)

}