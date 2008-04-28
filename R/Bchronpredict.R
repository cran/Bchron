Bchronpredict <- function(Bchrondata=list(SHOULDRUN=FALSE)) {

cat("Predict ages for the entire core. \n")
if(Bchrondata$SHOULDRUN == FALSE) {
    cat("No data found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

Bchrondata$chronsfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"chrons.txt", sep = "")
if(file.exists(Bchrondata$chronsfile)) {
  cat("Prediction stage appears to have been already run for this core. \n")
  cat("Do you wish to re-run? (y/n) \n")
  rerun <- scan(what = "", nlines = 1, quiet = TRUE)
  while(length(rerun)==0) rerun <- scan(what = "", nlines = 1, quiet = TRUE)
  if(rerun=="n" || rerun=="no") return(Bchrondata)
}

# Get the number of determinations
Temp <- read.table(Bchrondata$inputfile,header = TRUE)
Bchrondata$ndet <- nrow(Temp)

# Get the number of ddepths
Temp2 <- read.table(Bchrondata$ddepthfile, header = FALSE)
Bchrondata$nddepths <- nrow(Temp2)

# Get some useful files
Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "pars.txt",sep = "")
Bchrondata$outlierfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"outliers.txt", sep = "")

cat("Input date of core extraction in k cal yrs BP (leave blank for default of year 2000 = ",
    (1950 -2000)/1000,
    " k cal years BP):\n", sep = "")
Bchrondata$extractdate <- scan(what = "", nlines = 1, quiet = TRUE)
if (length(Bchrondata$extractdate) == 0) Bchrondata$extractdate <- (1950 - 2000)/1000
Bchrondata$numchron <- 10000

cat("Parameters file at", Bchrondata$parsfile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Output file at", Bchrondata$chronsfile, " \n")
cat("Design depths file at", Bchrondata$ddepthfile, " \n")
cat("Outlier output file at", Bchrondata$outlierfile, " \n")
cat("Number of determinations is", Bchrondata$ndet, " \n")
cat("Number of design depths is", Bchrondata$nddepths, " \n")
cat("Number of desired chronologies is", Bchrondata$numchron, " \n")
cat("Date of extraction of core (in k yrs BP) is", Bchrondata$extractdate,
    " \n")
cat("\n")
out <- .C("predict", as.character(Bchrondata$parsfile),
    as.character(Bchrondata$inputfile),
    as.character(Bchrondata$chronsfile),
    as.integer(Bchrondata$ndet),
    as.character(Bchrondata$ddepthfile),
    as.integer(Bchrondata$nddepths),
    as.integer(Bchrondata$numchron),
    as.double(Bchrondata$extractdate),
    as.character(Bchrondata$outlierfile))

cat("Completed!\n")
cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

Bchrondata$PREDICTEDALL <- TRUE
return(Bchrondata)

}