BchronreadEPD <- function() {

# Let's try and load in some html files direct from the EPD

cat("This function takes output taken directly from the EPD and converts \n")
cat("it for use in Bchron. For this to work, a file must be created which \n")
cat("contains only the dating information from the EPD. \n")

# Get the path
BADPATH <- TRUE
while(BADPATH==TRUE) {
  cat("Please enter the path to the directory that contains the EPD data file. \n")
  cat("Leave blank for the default of C:/Bchron/Input \n")
  EPDpath <- scan(what = "", nlines = 1, quiet = TRUE)
  if(length(EPDpath) == 0) EPDpath <- "C:/Bchron/Input"
  if(EPDpath==0) return(NULL)
  if(!file.exists(EPDpath)) {
    cat("Error: path cannot be found. \n")
    cat(paste("Currently it is set at: ",EPDpath,"\n", sep = ""))
    cat("Press <Enter> to try again...")
    readline()
    invisible()
  } else {
    BADPATH <- FALSE
  }
}

# Input file
BADINPUT <- TRUE
while(BADINPUT == TRUE) {
  cat("Now please enter the name of the EPD file you wish to use (including the file extension). \n")
  BADNAME <- TRUE
  while(BADNAME==TRUE) {
    EPDname <- scan(what = "", nlines = 1, quiet = TRUE)
    while(length(EPDname)==0) EPDname <- scan(what = "", nlines = 1, quiet = TRUE)
    if(length(EPDname)>0) BADNAME <- FALSE
    if(EPDname==0) return(NULL)
  }

  EPDfile <- paste(EPDpath,"/",EPDname,sep="")
  if(!file.exists(EPDfile)) {
      cat("Error: cannot find this file. \n")
      cat(paste("Check: ",EPDfile,"\n", sep = ""))
      cat("Press <Enter> to try again...")
      readline()
      invisible()
  } else {
    BADINPUT <- FALSE
  }
}

cat("Please enter a short name for the core. \n")
EPDfullname <- NULL
if(length(EPDfullname) == 0) EPDfullname <- scan(what = "", nlines = 1, quiet = TRUE,sep = "\t")
if(EPDfullname==0) return(NULL)
cat("Thank you.\n\n")

options(warn=-1)
EPDdata <- read.table(EPDfile,header=TRUE,sep="\t")
options(warn=1)

# Bchron input file should have:
# id	Age	sd	Depth	Thickness	Outlier1	Outlier2	Type

id <- EPDdata[,8]
age <- EPDdata[,5]
sd <- abs(EPDdata[,6]-EPDdata[,5])
depth <- EPDdata[,1]
thickness <- EPDdata[,2]
outlier1 <- rep(0.05,nrow(EPDdata))
outlier2 <- rep(0.001,nrow(EPDdata))
type <- rep(1,nrow(EPDdata))
BchronInput <- data.frame(id,age,sd,depth,thickness,outlier1,outlier2,type)

write.table(BchronInput,file=paste(EPDpath,"/",EPDfullname,".dat",sep=""),row.names=FALSE,sep="\t",quote=FALSE)

cat("Completed successfully. \n")
cat(paste("Bchron input created at: ",EPDpath,"/",EPDfullname,".dat \n",sep=""))
cat("Please now choose option 1 to use this data set. \n")
cat("Press <Enter> to continue...")
readline()
invisible()

}