BchronMCMC <-
function(Bchrondata,iterations=1e+05,burnin=10000,thinby=8,howmany=2000,defaults=FALSE,testrun=FALSE) {

cat("Setting up Bchron model run... \n")

if(file.exists(Bchrondata$parsfile)) {
    if(defaults!=TRUE) {
      cat("MCMC stage appears to have been already run for this core. \n")
      cat("Do you wish to re-run? (y/n) \n")
      rerun <- scan(what = "", nlines = 1, quiet = TRUE)
      while(length(rerun)==0) rerun <- scan(what = "", nlines = 1, quiet = TRUE)
      if(rerun=="n" || rerun=="no") return()
    }
}

# If a testrun, create an appropriate pars file
if(testrun==TRUE) {
    Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "parstestrun.txt",sep = "")
    iterations <- 100
    burnin <- 20
    howmany <- 10
    thinby <- 1
}

cat("Calibration curve at", Bchrondata$calibcurvefile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Output file at", Bchrondata$parsfile, " \n")
cat("Number of determinations is", nrow(Bchrondata$input), " \n")
cat("\n")
fails <- # Add a small amount of thickness (causes an error on linux)
thick <- Bchrondata$input[,5]
if(any(thick==0)) thick <- thick + 0.001

out <- .C("cpg", 
        as.double(Bchrondata$input[,2]/1000),
        as.double(Bchrondata$input[,3]/1000),
        as.double(Bchrondata$input[,4]/100),
        as.double(thick/100),
        as.double(Bchrondata$input[,6]),
        as.double(Bchrondata$input[,7]),
        as.integer(Bchrondata$input[,8]),
        as.character(Bchrondata$calibcurvefile),
        as.character(Bchrondata$parsfile),
        as.integer(nrow(Bchrondata$input)),
        as.integer(Bchrondata$bigcalsize),
        as.double(Bchrondata$lowcal),
        as.double(Bchrondata$highcal),
        as.integer(iterations),
        as.integer(burnin),
        as.integer(howmany),
        as.integer(thinby),
        as.integer(fails))

# See if the run was a success
Bchrondata$RUN <- as.logical(1-out[[12]][1])

if(Bchrondata$RUN==1) {

    if(testrun==FALSE) {
        cat("\n")
        cat("=======================\n")
        cat("=====Run completed=====\n")
        cat("=======================\n\n")
        cat(paste("Output completed: see ",Bchrondata$parsfile,sep=""),"\n")
    }

    if(testrun==TRUE) cat("Test run completed. Now do a longer run. \n")

} else {

    cat("\n")
    cat("======================\n")
    cat("======Run failed======\n")
    cat("======================\n")
    cat("\n")
        
}

}

