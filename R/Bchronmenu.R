`Bchronmenu` <-
function() {

Bchronversion <-"1.2"

cat("------------------------------- \n")
cat(paste("Welcome to Bchron version", Bchronversion, "\n"))
cat(paste("Author: Andrew Parnell, Trinity College Dublin\n"))
cat(paste("Please report bugs to: Andrew.Parnell@tcd.ie\n"))
cat("------------------------------- \n")

PATH <- NULL

EXIT <- FALSE
while(EXIT==FALSE)
{

choices <- c("Read in a Bchron data file (RUN THIS FIRST)","Calibrate a set of radiocarbon dates (30 seconds - 5 minutes)","Run Bchron (20 minutes - 12 hours)","Predict ages for the entire core (30 seconds - 5 minutes)","Predict ages for certain depths in the core (~1 minute)","Produce plots of the entire chronology (~2 mins)","Produce plots of the age for certain depths in the core (~30 seconds)","FIRST TIME USERS AND HELP SYSTEM","Exit")
title <- "The available options are:"
choose <- menu(choices,title = title)


#####################################################################################################

# Section 1
if(choose == 1) {

cat("You have chosen to read in a Bchron data file. \n")
cat("For this you will need to have created an appropriate file as outlined in the help section. \n")
cat("\n")

cat("Please enter the PATH to the file that contains the data file \n")
cat("eg C:\\cores\\mycore \n")
if(R.Version()$arch=="i386")
{
    cat("Leave blank for the default of C:\\Bchron \n")
} else {
    cat("Leave blank for the default of ~/Bchron \n")
}
cat("\n")
PATH <- scan(what="",nlines=1,quiet=TRUE)
if(length(PATH)==0 && R.Version()$arch=="i386") PATH <- "C:\\Bchron"
if(length(PATH)==0 && R.Version()$arch!="i386") PATH <- "~/Bchron"
cat("Now please enter the name of the .dat file you wish to run Bchron from. \n")
name <- scan(what="",nlines=1,quiet=TRUE)
cat("Please enter the full name of the core (leave blank if the same as above). \n")
fullname <- scan(what="",nlines=1,quiet=TRUE,sep="\t")
if(length(fullname)==0) fullname <- name
cat("Thank you.\n\n")

# CHECKS
if(R.Version()$arch=="i386")
{
    if(!file.exists(paste(PATH,"\\CalCurve\\BigCal.txt",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"\\CalCurve\\BigCal.txt\n",sep=""))
        stop("Calibration curve cannot be found",call.=FALSE)
    } 
    if(!file.exists(paste(PATH,"\\Input\\",name,".dat",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"\\Input\\",name,".dat\n",sep=""))
        stop("Data input file cannot be read",call.=FALSE)
    } 
} else {
    if(!file.exists(paste(PATH,"/CalCurve/BigCal.txt",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"/CalCurve/BigCal.txt\n",sep=""))
        stop("Calibration curve cannot be found",call.=FALSE)
    } 
    if(!file.exists(paste(PATH,"/Input/",name,".dat",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"/Input/",name,".dat\n",sep=""))
        stop("Data input file cannot be read",call.=FALSE)
    } 
}

cat("Press <Enter> to continue...")
readline()
invisible()
}

#####################################################################################################

if(choose == 2) {

cat("Now calibrating dates... \n")

if(length(PATH)==0) {
    cat("You have not entered any data. \n")
    cat("Please start again by running option 1. \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)

# Put in the path for the calibration curve
if(R.Version()$arch=="i386")
{
    calpath <- paste(PATH,"\\CalCurve\\BigCal.txt",sep="")
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"TrueDates.txt",sep="")
} else {
    calpath <- paste(PATH,"/CalCurve/BigCal.txt",sep="")
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    outfile <- paste(PATH,"/Output/",name,"TrueDates.txt",sep="")
}
calibrate(calpath,infile,outfile,ndet)


cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}

#####################################################################################################

if(choose == 3) {

cat("Setting up Bchron model... \n")

if(length(PATH)==0) {
    cat("You have not entered any data. \n")
    cat("Please start again by running option 1. \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)


# Put in the path for the calibration curve
if(R.Version()$arch=="i386")
{
    calpath <- paste(PATH,"\\CalCurve\\BigCal.txt",sep="")
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"pars.txt",sep="")
} else {
    calpath <- paste(PATH,"/CalCurve/BigCal.txt",sep="")
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    outfile <- paste(PATH,"/Output/",name,"pars.txt",sep="")
}
choices2 <- c("short","long","super long")

choose2 <- menu(choices2,title="Would you like a short or long run of the model?")

iterations <- 100000
burnin <- 10000
howmany <- 2000
thinby <- 8
if(choose2==2) {
    iterations <- 1000000
    burnin <- 200000
    howmany <- 20000
    thinby <- 75
}
if(choose2==3) {
    iterations <- 10000000
    burnin <- 2000000
    howmany <- 20000
    thinby <- 750
}


BchronMCMC(calpath,infile,outfile,ndet,iterations,burnin,howmany,thinby)

cat("\n Now checking convergence of Bchron parameters... \n")
pars <- read.table(outfile)
pars <-  pars[,c(seq(1,ndet),c(ncol(pars)-1,ncol(pars)))]

good <- try(boa.geweke(pars,0.1,0.5)[1,2],silent=TRUE)
if(is.numeric(good)) {
    pvals <- boa.geweke(pars,0.1,0.5)[,2]
    bad <- pvals<0.01
    vbad <- pvals<0.001

    if(sum(vbad[!is.na(vbad)])>0) {
    cat("SEVERE WARNING: it looks like",sum(vbad),"of the parameters have not converged. \n")
    } else if(sum(bad[!is.na(bad)])>0) {
    cat("WARNING:",sum(bad),"parameters may not have converged. \n")
    cat("Problem does not appear to be fatal. \n")
    }
    if(sum(bad[!is.na(bad)])>0) {
        cat("Plotting trace plots and densities. If these look unsatisfactory, do a longer run. \n")
        bad[is.na(bad)] <- FALSE
        badpars <- as.matrix(pars[,bad==TRUE])
        for(i in 1:ncol(badpars)) {
            if(i!=1) windows()
            par(mfrow=c(2,1))
            plot(badpars[,i],main="Trace plot",ylab=paste("Bad parameter",i))
            plot(density(badpars[,i]),main=paste("Density plot - bad parameter",i))
        }
        
        par(mfrow=c(1,1))

    }

} else {
    cat("Error in checking covergence of parameters. \n")
    cat("Check the CORENAMEpars.txt in the output directory for strange behaviour. \n")
}

cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}

#####################################################################################################

if(choose == 4) {

cat("Now predicting ages for the entire core \n")

if(length(PATH)==0) {
    cat("You have not entered any data. \n")
    cat("Please start again by running option 1. \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)

# Read in the ddepths file to get the number of ddepths
if(R.Version()$arch=="i386")
{
    Temp2 <- read.table(paste(PATH,"\\Input\\",name,"ddepths.txt",sep=""),header=FALSE)
} else {
    Temp2 <- read.table(paste(PATH,"/Input/",name,"ddepths.txt",sep=""),header=FALSE)
}
nddepths <- nrow(Temp2)

# Put in the paths for the various files
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    parsfile <- paste(PATH,"\\Output\\",name,"pars.txt",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"\\Input\\",name,"ddepths.txt",sep="")
    outlierfile <- paste(PATH,"\\Output\\",name,"Outliers.txt",sep="")
} else {
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    parsfile <- paste(PATH,"/Output/",name,"pars.txt",sep="")
    outfile <- paste(PATH,"/Output/",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"/Input/",name,"ddepths.txt",sep="")
    outlierfile <- paste(PATH,"/Output/",name,"Outliers.txt",sep="")
}
cat("Input date of core extraction in k cal yrs BP (leave blank for default = ",(1950-as.numeric(substr(date(),21,24)))/1000," k cal years BP):\n",sep="")
extract <- scan(what="",nlines=1,quiet=TRUE)
if(length(extract)==0) extract <- (1950-as.numeric(substr(date(),21,24)))/1000

numchron <- 10000

predictBchron(parsfile,infile,outfile,ndet,ddepthfile,nddepths,numchron,extract,outlierfile)

cat("Completed!\n")
cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()


}

#####################################################################################################

if(choose == 5) {

choices4 <- c("Fixed depths","Depth intervals")
title4 <- "Do you wish to use fixed depths in the core or would you rather use intervals?"
choose4 <- menu(choices4,title = title4)

if(length(PATH)==0) {
    cat("You have not entered any data. Make sure you choose option 1 with every run of Bchronmenu(). \n")
    cat("Please start again by running option 1. \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)

# Put in the paths for the various files
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    parsfile <- paste(PATH,"\\Output\\",name,"pars.txt",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"eventages.txt",sep="")
    outlierfile <- paste(PATH,"\\Output\\",name,"Outliers.txt",sep="")
} else {
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    parsfile <- paste(PATH,"/Output/",name,"pars.txt",sep="")
    outfile <- paste(PATH,"/Output/",name,"eventages.txt",sep="")
    outlierfile <- paste(PATH,"/Output/",name,"Outliers.txt",sep="")
}

if(!exists("extract")) {
   cat("Input date of core extraction in k cal yrs BP (leave blank for default = ",(1950-as.numeric(substr(date(),21,24)))/1000," k cal years BP):\n",sep="")
   extract <- scan(what="",nlines=1,quiet=TRUE)
   if(length(extract)==0) extract <- (1950-as.numeric(substr(date(),21,24)))/1000
}

####### Fixed depths
if(choose4==1) {

cat("Please create a file of the desired depths you wish to create \n")
cat("ages for in the input directory. \n")
cat("Make sure it is called MyCoreEventDepthsFixed.txt where 'MyCore' \n")
cat("is the name of the core. \n")
cat("Press <Enter> to continue...")
readline()
invisible()


# Read in the ddepths file to get the number of ddepths
if(R.Version()$arch=="i386")
{
    Temp2 <- read.table(paste(PATH,"\\Input\\",name,"EventDepthsFixed.txt",sep=""),header=FALSE)/100
} else {
    Temp2 <- read.table(paste(PATH,"/Input/",name,"EventDepthsFixed.txt",sep=""),header=FALSE)/100
}
nevents <- nrow(Temp2)

# Put in the paths for the various files
if(R.Version()$arch=="i386")
{
    evfile <- paste(PATH,"\\Input\\",name,"EventDepthsFixed.txt",sep="")
} else {
    evfile <- paste(PATH,"/Input/",name,"EventDepthsFixed.txt",sep="")
}

numchron <- 10000

predictBchron(parsfile,infile,outfile,ndet,evfile,nevents,numchron,extract,outlierfile)

cat("Completed!\n")
cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()


}


######### Random depths

if(choose4 == 2) {

cat("Please create a file of the desired depths interval you \n")
cat("wish to create ages for in the input directory. \n")
cat("The format of the file should look like the following: \n")
cat("261 265 \n")
cat("267 269 \n")
cat("271 282 \n")
cat("So each line is a single interval with the depths given in cm. \n")
cat("You can have as many intervals as you like for each event \n")
cat("Make sure it is called MyCoreEventDepthsEventName.txt where 'MyCore' \n")
cat("is the name of the core, and EventName is the name of the event in \n")
cat("question (eg GlendaloughEventDepthsUlmusDecline.txt). \n")
cat("Press <Enter> to continue...")
readline()
invisible()


cat("Now please enter the event name that you have given the file (eg UlmusDecline) \n")
badevname <- TRUE
while(badevname == TRUE) {
    evname <- scan(what="",nlines=1,quiet=TRUE)
    
    if(length(evname)>0) {
        # Read in the ddepths file to get the number of ddepths
        if(R.Version()$arch=="i386")
        {
            if(file.exists(paste(PATH,"\\Input\\",name,"EventDepths",evname,".txt",sep=""))) {
                Temp2 <- read.table(paste(PATH,"\\Input\\",name,"EventDepths",evname,".txt",sep=""),header=FALSE)/100
                badevname <- FALSE
            } else {
                cat(paste(PATH,"\\Input\\",name,"EventDepths",evname,".txt: file cannot be found. \n",sep=""))
            }
        } else {
            if(file.exists(paste(PATH,"/Input/",name,"EventDepths",evname,".txt",sep=""))) {
                Temp2 <- read.table(paste(PATH,"\\Input\\",name,"EventDepths",evname,".txt",sep=""),header=FALSE)/100
                badevname <- FALSE
            } else {
                cat(paste(PPATH,"/Input/",name,"EventDepths",evname,".txt: file cannot be found. \n",sep=""))
            }
        }
        nevents <- max(Temp2[,1])
    }
}


# Put in the paths for the various files
if(R.Version()$arch=="i386")
{
    outfile <- paste(PATH,"\\Output\\",name,"eventages",evname,".txt",sep="")
} else {
    outfile <- paste(PATH,"\\Output\\",name,"eventages",evname,".txt",sep="")
}

lowdepths <- Temp2[,1]
highdepths <- Temp2[,2]

numchron <- 10000
predictrandBchron(parsfile,infile,outfile,lowdepths,highdepths,length(lowdepths),ndet,numchron,extract,outlierfile)

cat("Completed!\n")
cat("\n \n \n")

cat("Press <Enter> to continue...")
readline()
invisible()

}

}

#####################################################################################################

if(choose == 6) {

cat("Now plotting chronology. \n")

if(length(PATH)==0) {
    cat("You have not entered any data. \n")
    cat("Please start again by running option 1. \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

# Create arguments
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Output\\",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"\\Input\\",name,"ddepths.txt",sep="")
    detsfile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    datesfile <- paste(PATH,"\\Output\\",name,"TrueDates.txt",sep="")
    rangesfile <- paste(PATH,"\\Output\\",name,"Ranges.txt",sep="")
} else {
    infile <- paste(PATH,"/Output/",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"/Input/",name,"ddepths.txt",sep="")
    detsfile <- paste(PATH,"/Input/",name,".dat",sep="")
    datesfile <- paste(PATH,"/Output/",name,"TrueDates.txt",sep="")
    rangesfile <- paste(PATH,"/Output/",name,"Ranges.txt",sep="")
}

choices3 <- c("colour","black and white")
choose3 <- menu(choices3,title="How would you like the chronology to be drawn?")

cols <- ifelse(choose3==1,TRUE,FALSE)

PlotBchron(fullname,infile,ddepthfile,detsfile,datesfile,rangesfile,cols,Bchronversion)
cat("Completed!\n")
cat("\n \n \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}

#####################################################################################################

if(choose == 7) {

cat("To use this option, you must have previously run option 5 \n")
cat("\n")
if(!exists("choose4")) {
    cat("Option 5 does not appear to have been run. \n")
    cat("Please start again by running option 1 and then running option 5. \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

cat(paste("Now plotting date estimates for each of the depths in ",name,"eventages.txt \n",sep=""))

if(length(PATH)==0) {
    cat("You have not entered any data. Make sure you choose option 1 with every run of Bchronmenu() \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
    Bchronmenu()
}

if(choose4==2) {
    badnum <- -1
    while(badnum<0) {
        cat("Please choose the event number you wish to plot \n")
        evnum <- as.integer(scan(what="",nlines=1,quiet=TRUE))       
        if(evnum<=nevents) {
            badnum <- 1
        } else {
            cat("Number invalid, enter a number between 1 and",nevents,"\n")
        }
    }
}    

# Create arguments
if(choose4==1) {
    if(R.Version()$arch=="i386")
    {
        infile <- paste(PATH,"\\Output\\",name,"eventages.txt",sep="")
        depthfile <- paste(PATH,"\\Input\\",name,"EventDepths.txt",sep="")
        outfile <- paste(PATH,"\\Output\\",name,"EventHDRs.txt",sep="")
    } else {
        infile <- paste(PATH,"/Output/",name,"eventages.txt",sep="")
        depthfile <- paste(PATH,"/Input/",name,"EventDepths.txt",sep="")
        outfile <- paste(PATH,"/Output/",name,"EventHDRs.txt",sep="")
    }
    PlotDens(infile,depthfile,fullname,outfile,Bchronversion,ngroup=0)
}
if(choose4==2) {
    if(R.Version()$arch=="i386")
    {
        infile <- paste(PATH,"\\Output\\",name,"eventages",evnum,".txt",sep="")
        depthfile <- paste(PATH,"\\Input\\",name,"EventDepthsIntervals.txt",sep="")
        outfile <- paste(PATH,"\\Output\\",name,"EventHDRs",evnum,".txt",sep="")
    } else {
        infile <- paste(PATH,"/Output/",name,"eventages",evnum,".txt",sep="")
        depthfile <- paste(PATH,"/Input/",name,"EventDepthsIntervals.txt",sep="")
        outfile <- paste(PATH,"/Output/",name,"EventHDRs",evnum,".txt",sep="")
    }
    PlotDens(infile,depthfile,fullname,outfile,Bchronversion,ngroup=evnum)
}

cat("\n")
cat("Done! \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}

#####################################################################################################

if(choose == 8) {

cat("Thank you for downloading Bchron, a (reasonably) user-friendly program for producing \n")
cat("reliably dated radiocarbon depth chronologies in Windows. If you are ever stuck whilst using this \n")
cat("program, you can always type help(Bchronmenu) at the command prompt for more information.\n")
cat("\n")
cat("WARNING: Some of these tasks, especially a long run of the Bchron model, may take several hours \n")
cat("Be aware before you start a long run that this may occur. It is often sensible to run these \n")
cat("overnight. The good thing is you only have to do it once for each core. \n")
cat("\n")
cat("To start, it is recommended for basic users to create a folder called Bchron at the \n")
cat("root directory of the hard disk. Within, you should create three directories, called inputs, \n")
cat("outputs and calcurve. The inputs directory should contain all the information about the cores for \n")
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
cat("3. (optional) 'coreeventdepths.txt', a single column file containing depths at which \n")
cat("   events of specific interest may have occurred. Bchron allows the ages of these to \n")
cat("   be looked at individually.\n")
cat("\n")
cat("Finally, the supplied 'BigCal.txt' file should be placed in the CalCurve directory. \n")
cat("Note that Bchron only uses (at present) the Northern hemisphere 2004 Intcal calibration \n")
cat("curve. Other calibration curves are not supported.\n")
cat("\n")
cat("With all this setup, it should be possible to follow the instructions after typing Bchronmenu() \n")
cat("at the command prompt. If you find any bugs, or wish to suggest enhancements, please contact \n")
cat("the author at Andrew.Parnell@tcd.ie \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}

if(choose == 9)
{
cat("Thank you. Exiting... \n")
EXIT=TRUE
}


}

}

