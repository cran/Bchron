Bchroncheck <-
function(Bchrondata) {

cat("Checking files for possible errors... \n\n")

if(!file.exists(Bchrondata$inputfile)) {
    cat("No data found. \n")
    cat("Please start again by running option 1, or calling Bchronload(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
}

# Checks:
# 1) Check that all files exist
# 2) Check 14C ages don't go outside the range of the calibration curves
# 3) Check input files make sense
# 4) Check ddepth file makes sense

# 1) Check files exist - calibcurvefile, inputfile, ddepthsfile
cat("1. Checking files exist...\n")
cat("a. Calibration curve......")
if(file.exists(Bchrondata$calibcurvefile)) {
  cat("done\n")
} else {
  cat("error\n")
  cat(paste("Calibration curve not found:",Bchrondata$calibcurvefile,"\n"))
}
cat("b. Determinations file......")
if(file.exists(Bchrondata$inputfile)) {
  cat("done\n")
} else {
  cat("error\n")
  cat(paste("Determinations file not found:",Bchrondata$inputfile,"\n"))
  cat("Press <Enter> to continue...")
  readline()
  invisible()
}

cat("\n")
# Check input file for correct format
cat("2. Checking determinations file for correct format...\n")
InputTemp <- read.table(Bchrondata$inputfile,header=TRUE)
cat("a. Number of columns.....")
if(ncol(InputTemp)==8) {
  cat("done\n")
} else {
  cat("error\n")
  cat(paste("Input file has incorrect number of columns:",Bchrondata$inputfile,"\n"))
  cat("Press <Enter> to continue...")
  readline()
  invisible()
}
cat("b. Checking for impossible depths.....")
tempdepths <- InputTemp[,4]
tempthick <- InputTemp[,5]
if(any(tempthick[duplicated(tempdepths)]==0)) {
  cat("error\n")
  cat(paste("Depths appear to be duplicated with zero thickness:",Bchrondata$inputfile,"\n"))
  cat("Press <Enter> to continue...")
  readline()
  invisible()
} else {
  cat("done\n")
}
cat("c. Checking outlier probabilities.....")
tempout <- c(InputTemp[,6],InputTemp[,7])
if(any(tempout>1) || any(tempout<0)) {
  cat("error\n")
  cat(paste("Outliers outside (0,1) range:",Bchrondata$inputfile,"\n"))
  cat("Press <Enter> to continue...")
  readline()
  invisible()
} else {
  cat("done\n")
}
cat("d. Checking for non-numeracy.....")
tempnums <- c(InputTemp[,2],InputTemp[,3],InputTemp[,4],InputTemp[,5],InputTemp[,6],InputTemp[,7],InputTemp[,8])
if(!is.numeric(tempnums)) {
  cat("error\n")
  cat(paste("Check input file for non-numeric values:",Bchrondata$inputfile,"\n"))
  cat("Press <Enter> to continue...")
  readline()
  invisible()
} else {
  cat("done\n")
}

# 3) Check 14C ages don't go outside the range of the calibration curve
cat("\n")
cat("3. Checking calibration curve.....")
CalCurveTemp <- read.table(Bchrondata$calibcurvefile,header=FALSE)
Dates <- InputTemp[,2]
Sds <- InputTemp[,3]
Types <- InputTemp[,8]
C14dates <- Dates[Types==1]
if(length(C14dates)==0) {
    cat("\n Notice: no 14C dates found - proceed without calibration curve checking. \n")
} else {
    C14sds <- Sds[Types==1]
    lowdates <- round(C14dates-3*C14sds,0)
    highdates <- round(C14dates+3*C14sds,0)
    if(any((lowdates-5)<0) || any((highdates-5)>length(CalCurveTemp[,1])) ) {
      cat("error\n")
      cat("14C dates are outside the range of the calibration curve. \n")
      cat(paste("Check:",Bchrondata$inputfile,"and",Bchrondata$calibcurvefile,"\n"))
      cat("Press <Enter> to continue...")
      readline()
      invisible()
   } else {

      lowcalagelookup <- CalCurveTemp[lowdates-5,1]
      highcalagelookup <- CalCurveTemp[highdates-5,1]
      if(any(is.na(lowcalagelookup)) || any(is.na(highcalagelookup)) ) {
         cat("error\n")
         cat("14C dates are outside the range of the calibration curve. \n")
         cat(paste("Check:",Bchrondata$inputfile,"and",Bchrondata$calibcurvefile,"\n"))
         cat("Press <Enter> to continue...")
         readline()
         invisible()
      } else {
         cat("done\n")
      }
  }

}

# Check output files for mis-matches (if they exist)
if(file.exists(Bchrondata$chronsfile)) {
    cat("\n")
    cat("4. Checking output files.....")
    chrons <- read.table(Bchrondata$chronsfile)
    if(ncol(chrons)!=length(Bchrondata$outdepths)) {
      cat("warning\n")
      cat("Number of output depths does not match size of chronology. \n")
      cat("Re-run prediction stage. \n")
      cat("Press <Enter> to continue...")
      readline()
      invisible()
    } else {
      cat("done\n")
}
    
    
    

}


cat("\n")
cat("=====================\n")
cat("All checks completed. \n")
cat("=====================\n")
cat("Press <Enter> to continue...")
readline()
invisible()

}
