Bchroncalibrate <-
function(Bchrondata,iterations=500000,burnin=50000,thinby=45,howmany=50000,defaults=FALSE) {

cat("Calibrating radiocarbon dates... \n")

# Create output file path
if(file.exists(Bchrondata$calibdatesfile)) {
    if(defaults!=TRUE) {
      cat("Calibration stage appears to have been already run for this core. \n")
      cat("Do you wish to re-run? (y/n) \n")
      rerun <- scan(what = "", nlines = 1, quiet = TRUE)
      while(length(rerun)==0) rerun <- scan(what = "", nlines = 1, quiet = TRUE)
      if(rerun=="n" || rerun=="no") return()
    }
}

cat("Calibration curve at", Bchrondata$calibcurvefile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Output file at", Bchrondata$calibdatesfile, " \n")
cat("Number of determinations is", nrow(Bchrondata$input), " \n")
cat("\n")

# Give calibrate code file c14dates and c14errors only

out <- .C("calibrate", 
        as.double(Bchrondata$input[,2]/1000),
        as.double(Bchrondata$input[,3]/1000),
        as.integer(Bchrondata$input[,8]),
        as.character(Bchrondata$calibcurvefile),
        as.character(Bchrondata$calibdatesfile),
        as.integer(nrow(Bchrondata$input)),
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

# Create a calibrated dates ranges file
Bchrondata$calibdates <- read.table(Bchrondata$calibdatesfile)
calibranges <- matrix(0,ncol=3,nrow=ncol(Bchrondata$calibdates))
cat(paste("Date","2.5%","50%","97.5%","\n"),file=Bchrondata$calibrangesfile,append=FALSE)
for(i in 1:ncol(Bchrondata$calibdates)) {
    calibranges[i,] <- quantile(Bchrondata$calibdates[,i],probs=c(0.025,0.5,0.975))
    cat(c(paste(Bchrondata$input[i,1]),round(calibranges[i,],3),"\n"),file=Bchrondata$calibrangesfile,append=TRUE)
}

}
