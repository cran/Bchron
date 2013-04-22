Bchronmenu <-
function ()
{
EXIT <- FALSE

Bchronver <- read.dcf(file=system.file("DESCRIPTION", package="Bchron"),fields="Version")

cat("------------------------------- \n")
cat(paste("Welcome to Bchron version",Bchronver,"\n"))
cat("Author: Andrew Parnell, University College Dublin\n")
cat("Please report bugs to: Andrew.Parnell@ucd.ie\n")
cat("------------------------------- \n \n")
cat("Option 1 must be chosen with every new call to Bchronmenu() \n")
cat("Option 0 will return you to the main menu, or \n")
cat("press <Esc> to exit the program at any time. \n \n")

while(EXIT==FALSE)
{

choices <- c("Load in or set up a Bchron data file",
      "Calibrate a set of radiocarbon dates",
      "Run Bchron",
      "Predict and plot ages for the entire core",
      "Predict and plot ages for an event in the core",
      "Exit")
title <- "The available options are:"
choose <- menu(choices, title = title)

# Load in data
if(choose == 1) {
  
    path <-  dlgDir(title="Please select your Bchron directory")$res    
  
#    BADPATH <- TRUE
#    while(BADPATH==TRUE) {
#      cat("Please enter the path to the Bchron directory that contains the data files \n")
#      cat("and the calibration curve. \n")
#      cat("Leave blank for the default of",getwd(),"\n")
#      path <- scan(what = "", nlines = 1, quiet = TRUE)
#      if(length(path) == 0) path <- getwd()
#      if(path==0) break
#      if(!file.exists(path)) {
#        cat("Error: path cannot be found. \n")
#        cat(paste("Currently it is set at: ",path,"\n", sep = ""))
#        cat("Press <Enter> to try again...")
#        readline()
#        invisible()
#      } else {
#        BADPATH <- FALSE
#      }
#    }
#    if(path==0) next

    cat("Please select an available (.bch) calibration curve: \n")
    choices2 <- c(list.files(paste(path,"/CalCurve/",sep=""),pattern="bch$"))
    title2 <- "List of available calibration curves:"
    choose2 <- menu(choices2, title = title2)
    if(choose2==0) next
    calibname <- strsplit(choices2[choose2],".",fixed=TRUE)[[1]][1]

    cat("Please select an available input (.dat) file \n")
    choices3 <- c(list.files(paste(path,"/Input/",sep=""),pattern="dat$"))
    title3 <- "List of available input files:"
    choose3 <- menu(choices3, title = title3)
    if(choose3==0) next
    inputfile <- paste(path,"/Input/",choices3[choose3],sep="")
    name <- strsplit(choices3[choose3],".",fixed=TRUE)[[1]][1]
    
    cat("Please enter the top depth (in cm) at which chronology age is required \n")
    cat("(or leave blank for default of top dated depth to bottom dated depth). \n")
    topoutdepth <- scan(nlines = 1, quiet = TRUE)
    if(length(topoutdepth) == 0) {
        tmp <- read.table(inputfile,header=TRUE)
        outdepths <- seq(min(tmp[,4]),max(tmp[,4]),length=200)
    }
    
    if(length(topoutdepth)>0) {
        cat("Please enter the bottom depth (in cm) at which chronology age is required. \n")
        botoutdepth <- scan(nlines = 1, quiet = TRUE)
        while(length(topoutdepth) == 0) {
            botoutdepth <- scan(nlines = 1, quiet = TRUE)
        }
    
        cat("Please enter the step size between output depths. \n")
        stepdepth <- scan(nlines = 1, quiet = TRUE)
        while(length(stepdepth) == 0) {
            stepdepth <- scan(nlines = 1, quiet = TRUE)
        }
        outdepths <- seq(topoutdepth,botoutdepth,by=stepdepth)
    }
        
    # Get extraction date
    cat("Please input date of core extraction in k cal yrs BP \n")
    cat("(leave blank for default of year 2000 = ",(1950 -2000)/1000," k cal years BP):\n", sep = "")
    extractdate <- scan(what = "", nlines = 1, quiet = TRUE)
    if (length(extractdate) == 0) extractdate <- (1950 - 2000)/1000
    
    # Get fullname
    cat("Please enter the full name of the core (used for graph titles; leave blank if the same as above). \n")
    fullname <- scan(what = "", nlines = 1, quiet = TRUE,sep = "\t")
    if(length(fullname) == 0) fullname <- name
    if(fullname==0) next
    
    Bchronrun <- Bchronload(name,fullname,path,outdepths,calibname,extractdate,check=TRUE)   
}  
    

# Calibrate the 14C dates
if(choose==2) {
  
    cat("============================================================================\n")
    cat("Calibrate radiocarbon dates \n")
    cat("============================================================================\n")
    
    choices2 <- c("standard", "long", "super-long")
    choose2 <- menu(choices2, title = "What type of Bchron model run would you like?")  
    if(choose2==0) next
    if(choose2==1) Bchroncalibrate(Bchronrun)
    if(choose2==2) Bchroncalibrate(Bchronrun,iterations=1e+06,burnin=200000,thinby=80,howmany=50000)
    if(choose2==3) Bchroncalibrate(Bchronrun,iterations=2e+06,burnin=400000,thinby=160,howmany=50000)
    cat("Press <Enter> to continue or <Esc> to exit...")
    readline()
    invisible()
    
    cat("Do you wish to check convergence? (y/n) \n")
    checkconverge <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(checkconverge)==0) checkconverge <- scan(what = "", nlines = 1, quiet = TRUE)
    if(checkconverge==0) next
    if(checkconverge=="y" || checkconverge=="yes") {
      Bchronconvergecheck(Bchronrun,dates=TRUE)
    }
    cat("Do you wish to plot the calibrated dates? (y/n) \n")
    plotdates <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(plotdates)==0) plotdates <- scan(what = "", nlines = 1, quiet = TRUE)
    if(plotdates==0) next
    if(plotdates=="y" || plotdates=="yes") {
       Bchronplot(Bchronrun,dates=TRUE,chrons=FALSE)
    }
}

# Run the Bchron model
if(choose==3) {

    cat("============================================================================\n")
    cat("Run Bchron model \n")
    cat("WARNING: longer runs may take up to 12 hours to complete. \n")
    cat("============================================================================\n")
    cat("Press <Enter> to continue or <Esc> to exit...")
    readline()
    invisible()

    choices2 <- c("test run","standard", "long", "super-long")
    choose2 <- menu(choices2, title = "What type of Bchron model run would you like?")
    if(choose2==0) next  
    if(choose2==1) BchronMCMC(Bchronrun,testrun=TRUE)
    if(choose2==2) BchronMCMC(Bchronrun)
    if(choose2==3) BchronMCMC(Bchronrun,iterations=2e+05,burnin=20000,thinby=75,howmany=20000)
    if(choose2==4) BchronMCMC(Bchronrun,iterations=1e+07,burnin=2e+06,thinby=750,howmany=20000)
    cat("Press <Enter> to continue or <Esc> to exit...")
    readline()
    invisible()
    
    if(choose2>1) {
        cat("Do you wish to check convergence? (y/n) \n")
        checkconverge <- scan(what = "", nlines = 1, quiet = TRUE)
        while(length(checkconverge)==0) checkconverge <- scan(what = "", nlines = 1, quiet = TRUE)
        if(checkconverge==0) next
        if(checkconverge=="y" || checkconverge=="yes") {
            Bchronconvergecheck(Bchronrun)
        }
    }
}

# Predict for all depths in the core
if(choose==4) {
    Bchronpredict(Bchronrun)

    cat("Do you wish to plot the chronology? (y/n) \n")
    plotchron <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(plotchron)==0) plotchron <- scan(what = "", nlines = 1, quiet = TRUE)
    if(plotchron==0) next
    if(plotchron=="y" || plotchron=="yes") {
      Bchrondata <- Bchronplot(Bchronrun)
    }
    cat("Press <Enter> to continue or <Esc> to exit...")
    readline()
    invisible()
}

if(choose==5) {
    # Predict just an event
    cat("============================================================================\n")
    cat("For each event, please create a file of the desired depths intervals \n")
    cat("for which you wish to create ages. Put the files in the input directory. \n")
    cat("The format of the file should look like the following: \n")
    cat("261 265 \n")
    cat("267 269 \n")
    cat("271 282 \n")
    cat("So each line is a single interval with the depths given in cm. \n")
    cat("You can have as many intervals as you like for each event \n")
    cat("Make sure it is called MyCoreEventDepthsEventName.txt where 'MyCore' \n")
    cat("is the name of the core, and EventName is the name of the event in \n")
    cat("question (eg GlendaloughEventDepthsUlmusDecline.txt). \n")
    cat("============================================================================\n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    cat("Now please enter the event name. \n")
    event <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(event)==0) event <- scan(what = "", nlines = 1, quiet = TRUE)
    if(event==0) next
    Bchronpredictevent(Bchronrun,event=event) 
    cat("Press <Enter> to continue...")
    readline()
    invisible()
  
    cat("Do you wish to plot this event? (y/n) \n")
    plotevent <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(plotevent)==0) plotevent <- scan(what = "", nlines = 1, quiet = TRUE)
    if(plotevent==0) next
    if(plotevent=="y" || plotevent=="yes") {
        cat("Please enter the full name of the event in inverted commas (used for plotting - leave blank for default). \n")
        eventname <- scan(what = "", nlines = 1, quiet = TRUE)
        if(eventname==0) next
        if(length(eventname)==0) {
            eventname <- event
        } else {
            eventname <- as.character(eventname)
        }
        Bchronplotevent(Bchronrun,event=event,eventname=eventname)
    }
}

# Exit
if(choose==6) {
  cat("Thank you. Exiting... \n")
  EXIT <- TRUE
}

}

}
