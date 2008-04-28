Bchronconvergecheck <- function(Bchrondata=list(RUN=FALSE),CALDATES=FALSE) {

if (Bchrondata$RUN == FALSE) {
    cat("No data run found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

if(CALDATES==TRUE && Bchrondata$CALIBRATED==FALSE) {
    cat("No calibrated data run found. \n")
    cat("Please start again by running option 1, or calling Bchronloaddata(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    return(Bchrondata)
}

library(coda)

# Get the number of determinations
Temp <- read.table(Bchrondata$inputfile,header = TRUE)
Bchrondata$ndet <- nrow(Temp)

cat("\n")
if(CALDATES==FALSE) {
  cat("Checking Bchron model run convergence... \n")
  pars <- read.table(Bchrondata$parsfile)
  pars <- pars[,c(seq(1,Bchrondata$ndet),c(ncol(pars)-1,ncol(pars)))]
}
if(CALDATES==TRUE) {
  cat("Checking Bchron calibrated dates convergence... \n")
  pars <- read.table(Bchrondata$calibdatesfile)
}

if(nrow(pars)<2000) cat("\nWARNING: This is a very short run - convergence results may be meaningless.\n\n")


good <- try(geweke.diag(pars)[[1]],silent=TRUE)
if (is.numeric(good)) {
    cat("========================================================================================\n")
    cat("If many of these values are small (less than 0.01), a longer run will be required.      \n")
    cat("If any of these values are NA, there may have been a problem. The core should be re-run.\n")
    cat("========================================================================================\n\n")
    cat("Worst parameters are ... \n")
    temp <- geweke.diag(pars)[[1]]
    pvals <- c(pnorm(temp[temp<0]),1-pnorm(temp[temp>0]))
    print(sort(c(pnorm(temp[temp<0]),1-pnorm(temp[temp>0])))[1:min(10,ncol(pars))])
    bad <- pvals < 0.01
    vbad <- pvals < 0.001
    if (sum(vbad[!is.na(vbad)]) > 0) {
      cat("SEVERE WARNING: it looks like some of the parameters have not converged. Try a longer run. \n")
    } else if (sum(bad[!is.na(bad)]) > 0) {
      cat("WARNING:", sum(bad), "parameters may not have converged. \n")
      cat("Problem does not appear to be fatal. Try a longer run. \n")
    } else {
      cat("Results are satisfactory. \n")
    }
} else {
    cat("Error in checking covergence of parameters. \n")
    if(CALDATES==TRUE) cat("Check the file", Bchrondata$parsfile,"in the output directory for strange behaviour. \n")
    if(CALDATES==FALSE) cat("Check the file", Bchrondata$caldatesfile,"in the output directory for strange behaviour. \n")
    cat("Or try a longer run.\n")
}
cat("\n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}