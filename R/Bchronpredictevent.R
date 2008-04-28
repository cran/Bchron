Bchronpredictevent <- function(Bchrondata) {

cat("Predict ages for an event in the core. \n")
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

# Get some useful files
Bchrondata$parsfile <- paste(Bchrondata$path, "/Output/", Bchrondata$name, "pars.txt",sep = "")
Bchrondata$outlierfile <- paste(Bchrondata$path, "/Output/",Bchrondata$name,"outliers.txt", sep = "")

cat("Input date of core extraction in k cal yrs BP (leave blank for default of year 2000 = ",
    (1950 - 2000)/1000,
    " k cal years BP):\n", sep = "")
Bchrondata$extractdate <- scan(what = "", nlines = 1, quiet = TRUE)
if (length(Bchrondata$extractdate) == 0) Bchrondata$extractdate <- (1950 - 2000)/1000

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
cat("Press <Enter> to continue...")
readline()
invisible()

numevents <- 0
while(numevents < 1) {
  cat("How many events do you wish to input? \n")
  numevents <- scan(what = "integer", nlines = 1, quiet = TRUE)
  while(length(numevents)==0) numevents <- scan(what = "integer", nlines = 1, quiet = TRUE)
}
Bchrondata$eventnames <- rep("",times=numevents)
Bchrondata$eventfullnames <- rep("",times=numevents)
Bchrondata$eventdepthfile <- rep("",times=numevents)
Bchrondata$eventagefile <- rep("",times=numevents)

for(i in 1:numevents) {
  if(numevents>1) cat(paste("Event",i,"\n"))
  cat("Now please enter the event name that you have given the file (eg UlmusDecline) \n")
  badevname <- TRUE
  while (badevname == TRUE) {
    evname <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(evname)==0) evname <- scan(what = "", nlines = 1, quiet = TRUE)
    Bchrondata$eventdepthfile[i] <- paste(Bchrondata$path,"/Input/",Bchrondata$name,"EventDepths",evname,".txt",sep="")
    Bchrondata$eventagefile[i] <- paste(Bchrondata$path,"/Output/",Bchrondata$name,"EventAges",evname,".txt",sep="")
    cat("Please enter the full name of the event (eg Ulmus Decline; for plotting purposes) \n")
    evfullname <- scan(what = "character", nlines = 1, quiet = TRUE,sep=",")
    if(length(evfullname)==0) evfullname <- evname
    if(file.exists(Bchrondata$eventdepthfile[i])) {
      badevname <- FALSE
      Temp2 <- read.table(Bchrondata$eventdepthfile[i])/100
      Bchrondata$eventnames[i] <- evname
      Bchrondata$eventfullnames[i] <- evfullname
    } else {
      cat(paste(Bchrondata$path, "/Input/", Bchrondata$name, "EventDepths",evname, ".txt \n"),sep="")
      cat("File cannot be found. \n")
    }
  }
  
  lowdepths <- Temp2[, 1]
  highdepths <- Temp2[, 2]
  nddepthints <- length(lowdepths)
  numchron <- 10000
  
  cat("Parameters file at", Bchrondata$parsfile, " \n")
  cat("Input file at", Bchrondata$inputfile, " \n")
  cat("Event depths file at", Bchrondata$eventdepthfile[i], " \n")
  cat("Event ages file at", Bchrondata$eventagefile[i], " \n")
  cat("Range of design depths is \n")
  for(j in 1:nddepthints) {
      cat(lowdepths[j],"to",highdepths[j],"\n")
  }
  cat("Outlier output file at", Bchrondata$outlierfile, " \n")
  cat("Number of determinations is", Bchrondata$ndet, " \n")
  cat("Number of design depths intervals is", nddepthints, " \n")
  cat("Number of desired chronologies is", numchron, " \n")
  cat("Date of extraction of core (in k yrs BP) is", Bchrondata$extractdate, " \n")

  cat("\n")

  out <- .C("predictrand",
    as.character(Bchrondata$parsfile),
    as.character(Bchrondata$inputfile),
    as.character(Bchrondata$eventagefile[i]),
    as.double(lowdepths),
    as.double(highdepths),
    as.integer(nddepthints),
    as.integer(Bchrondata$ndet),
    as.integer(numchron),
    as.double(Bchrondata$extractdate),
    as.character(Bchrondata$outlierfile))

  cat("Completed!\n")
  cat("\n \n \n")
  cat("Press <Enter> to continue...")
  readline()
  invisible()
  
}

  return(Bchrondata)
}