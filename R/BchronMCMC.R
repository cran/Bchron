BchronMCMC <- function(Bchrondata,iterations=0,burnin=0,thinby=1,howmany=1) {

# Get coda
library(coda)

cat("Setting up Bchron model... \n")
if (Bchrondata$SHOULDRUN == FALSE) {
    cat("No data found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

# Get the number of determinations
Temp <- read.table(Bchrondata$inputfile,header = TRUE)
Bchrondata$ndet <- nrow(Temp)

# Get an appropriate output file
Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "pars.txt",sep = "")

choose2 <- 2
if(iterations==0) {

# Determine the run size
cat("====================================================================\n")
cat("WARNING: running the model here will over-write previous model runs. \n")
cat("WARNING: longer runs may take up to 12 hours to complete. \n")
cat("====================================================================\n")
cat("Select 0 to exit this menu. \n")
choices2 <- c("test run","short ", "long", "super-long")
choose2 <- menu(choices2, title = "What type of Bchron model run would you like?")
iterations <- 100
burnin <- 20
howmany <- 10
thinby <- 1
if(choose2==1) Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "parstestrun.txt",sep = "")
if(choose2 == 2) {
    iterations <- 1e+05
    burnin <- 10000
    howmany <- 2000
    thinby <- 8
}
if (choose2 == 3) {
    iterations <- 1e+06
    burnin <- 2e+05
    howmany <- 20000
    thinby <- 75
}
if (choose2 == 4) {
    iterations <- 1e+07
    burnin <- 2e+06
    howmany <- 20000
    thinby <- 750
}
if(choose2 == 0) return(Bchrondata)
}

cat("Calibration curve at", Bchrondata$calibcurvefile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Output file at", Bchrondata$parsfile, " \n")
cat("Number of determinations is", Bchrondata$ndet, " \n")
cat("\n")
fails <- 0
out <- .C("cpg", as.character(Bchrondata$calibcurvefile),
        as.character(Bchrondata$inputfile),
        as.character(Bchrondata$parsfile),
        as.integer(Bchrondata$ndet),
        as.integer(Bchrondata$bigcalsize),
        as.double(Bchrondata$lowcal),
        as.double(Bchrondata$highcal),
        as.integer(iterations),
        as.integer(burnin),
        as.integer(howmany),
        as.integer(thinby),
        as.integer(fails))

# Only check parameters if run was a success
Bchrondata$RUN <- as.logical(1-out[[12]][1])

if(Bchrondata$RUN==1) {

if(choose2 > 1) {
cat("\n")
cat("=======================\n")
cat("=====Run completed=====\n")
cat("=======================\n\n")
cat(paste("Output completed: see ",Bchrondata$parsfile,sep=""))
}

if(choose2 == 1) cat("Test run completed. Now do a longer run.")

cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

return(Bchrondata)

} else {
  cat("\n")
  cat("======================\n")
  cat("======Run failed======\n")
  cat("======================\n")
  cat("\n")
  cat("Press <Enter> to continue...")
  readline()
  invisible()
  return(Bchrondata)
}

}