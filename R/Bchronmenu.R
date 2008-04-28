Bchronmenu <- function ()
{

Bchrondata <- list()
Bchrondata$SHOULDRUN <- FALSE
Bchrondata$EXIT <- FALSE
Bchrondata$CALIBRATED <- FALSE
Bchrondata$RUN <- FALSE
Bchrondata$version <- "2.0"

# Need this library
library(coda)
library(hdrcde)

cat("------------------------------- \n")
cat(paste("Welcome to Bchron version", Bchrondata$version, "\n"))
cat(paste("Author: Andrew Parnell, Trinity College Dublin\n"))
cat(paste("Please report bugs to: Andrew.Parnell@tcd.ie\n"))
cat("------------------------------- \n \n")
    
cat("Tip: press 0 at any time to return to the main menu, \n")
cat("or Esc to exit the program. \n \n")

while(Bchrondata$EXIT==FALSE)
{

choices <- c("Load in or set up a Bchron data file",
      "Calibrate a set of radiocarbon dates",
      "Run Bchron",
      "Predict and plot ages for the entire core",
      "Predict and plot ages for certain events in the core",
      "Change calibration curve used",
      "Check data files for possible errors",
      "First time users and help system",
      "Exit")
title <- "The available options are:"
choose <- menu(choices, title = title)

# Load in data
if (choose == 1) {
  Bchrondata <- Bchronloaddata(Bchrondata$version)
  if(Bchrondata$SHOULDRUN==TRUE) {
    cat("Do you wish to check the data are in the correct format? (y/n) \n")
    docheck <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(docheck)==0) docheck <- scan(what = "", nlines = 1, quiet = TRUE)
    if(docheck=="y" || docheck=="yes") {
      Bchrondata <- Bchroncheck(Bchrondata)
    }
  }
}
# Calibrate the 14C dates
if (choose == 2) Bchrondata <- Bchroncalibrate(Bchrondata)

# Run the Bchron model
if (choose == 3) {
  Bchrondata <- BchronMCMC(Bchrondata)
  if(Bchrondata$SHOULDRUN==TRUE) {
    cat("Do you wish to check convergence? (y/n) \n")
    checkconverge <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(checkconverge)==0) checkconverge <- scan(what = "", nlines = 1, quiet = TRUE)
    if(checkconverge=="y" || checkconverge=="yes") {
      Bchronconvergecheck(Bchrondata)
    }
  }
}

# Predict for all depths in the core
if (choose == 4) {
  Bchrondata <- Bchronpredict(Bchrondata)
  if(Bchrondata$SHOULDRUN==TRUE) {
    cat("Do you wish to plot the chronology? (y/n) \n")
    plotchron <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(plotchron)==0) plotchron <- scan(what = "", nlines = 1, quiet = TRUE)
    if(plotchron=="y" || plotchron=="yes") {
      Bchrondata <- Bchronplotter(Bchrondata)
    }
  }
}

if (choose == 5) {
  # Predict just an event
  Bchrondata <- Bchronpredictevent(Bchrondata)
  
  if(Bchrondata$SHOULDRUN==TRUE) {
    cat("Do you wish to plots these events? (y/n) \n")
    plotevent <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(plotevent)==0) plotevent <- scan(what = "", nlines = 1, quiet = TRUE)
    if(plotevent=="y" || plotevent=="yes") {
      Bchrondata <- Bchronplotevent(Bchrondata)
    }
  }
}

if (choose == 6) {
# Change calibration curve
   Bchrondata <- Bchronchangecalcurve(Bchrondata)
   if(Bchrondata$SHOULDRUN==TRUE) {
    cat("Do you wish to check the data agree with this calibration curve? (y/n) \n")
    docheck <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(docheck)==0) docheck <- scan(what = "", nlines = 1, quiet = TRUE)
    if(docheck=="y" || docheck=="yes") {
      Bchrondata <- Bchroncheck(Bchrondata)
    }
  }
}

# Run a data check
if(choose == 7) {
  Bchroncheck(Bchrondata)
}

# Some help
if(choose == 8) {
    cat("Thank you for downloading Bchron, a (reasonably) user-friendly program for producing \n")
    cat("reliably dated radiocarbon depth chronologies. If you are ever stuck whilst using this \n")
    cat("program, you can always type help(Bchronmenu) at the command prompt for more information.\n")
    cat("\n")
    cat("WARNING: Some of these tasks, especially a long run of the Bchron model, may take several hours. \n")
    cat("Be aware before you start a long run that this may occur. It is often sensible to run these \n")
    cat("overnight. The good thing is you only have to do it once for each core. \n")
    cat("\n")
    cat("To start, it is recommended for basic users to create a folder called Bchron at the \n")
    cat("root directory of the hard disk. Within, you should create three directories, called input, \n")
    cat("output and calcurve. The inputs directory should contain all the information about the cores for \n")
    cat("which you require chronologies. This should include:\n")
    cat("1. 'core.dat', a tab-delimited file with columns for the laboratory code of the sample, the \n")
    cat("   14C age, the error, the depth (in cm), the thickness (in cm), the probability of \n")
    cat("   being a standard outlier, the probability of being an extreme outlier and the type of date information. \n")
    cat("   Note that there are three allowed types: 1) a standard radiocarbon date, 2) a calendar age with a Normally \n")
    cat("   distributed error, 3) a calendar age with a Uniformly distributed error. For types 1) and 2), the error given \n")
    cat("   in column 3 is a 1-sigma standard error. For type 3), the error is the distance to the upper or lower bound \n")
    cat("   of the desired Uniform distibution. The format of the .dat file can \n")
    cat("   be copied from the supplied Glendalough.dat file. Note that the two outlier columns can \n")
    cat("   generally be left as is, though you may want to set some of them to 0 if the supplied dates are \n")
    cat("   not radiocarbon ages.\n")
    cat("2. 'coreddepths.txt', a single column file which contains each of the desired depths, for \n")
    cat("   example, where pollen was collected.\n")
    cat("3. (optional) 'coreeventdepthseventname.txt', a 2-column file containing the depth ranges \n")
    cat("   at which an event of specific interest may have occurred. \n")
    cat("\n")
    cat("Finally, the supplied 'IntCal04.bch' file should be placed in the CalCurve directory. \n")
    cat("Note that Bchron only uses (at present) the Northern hemisphere 2004 Intcal calibration \n")
    cat("curve. Other calibration curves can be downloaded from eg http://www.radiocarbon.org/IntCal04.htm. \n")
    cat("They then need to be converted via the Bchron menu system or by the function Bchronchangecalcurve().\n")
    cat("\n")
    cat("With all this setup, it should be possible to follow the instructions after typing Bchronmenu() \n")
    cat("at the command prompt. If you find any bugs, or wish to suggest enhancements, please contact \n")
    cat("the author at Andrew.Parnell@tcd.ie \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
}

# Exit
if (choose == 9) {
  cat("Thank you. Exiting... \n")
  Bchrondata$EXIT = TRUE
}

}

}