Bchronpredict <- function(Bchrondata,numchron=10000,defaults=FALSE) {

cat("Predict ages for the entire core. \n")

if(!file.exists(Bchrondata$parsfile)) {
    cat("No model run found. \n")
    cat("Please start again by running option 1, or calling Bchronload(). \n \n")
    cat("Press <Enter> to continue...")
    readline()
    invisible()
}

if(defaults==FALSE) {
    if(file.exists(Bchrondata$chronsfile)) {
        cat("Prediction stage appears to have been already run for this core. \n")
        cat("Do you wish to re-run? (y/n) \n")
        rerun <- scan(what = "", nlines = 1, quiet = TRUE)
        while(length(rerun)==0) rerun <- scan(what = "", nlines = 1, quiet = TRUE)
        if(rerun=="n" || rerun=="no") return()
    }
}

cat("Parameters file at", Bchrondata$parsfile, " \n")
cat("Input file at", Bchrondata$inputfile, " \n")
cat("Output file at", Bchrondata$chronsfile, " \n")
cat("Outlier output file at", Bchrondata$outlierfile, " \n")
cat("Number of determinations is", nrow(Bchrondata$input), " \n")
cat("Number of design depths is", length(Bchrondata$outdepths), " \n")
cat("Number of desired chronologies is", numchron, " \n")
cat("Date of extraction of core (in k yrs BP) is", Bchrondata$extractdate,
    " \n")
cat("\n")

out <- .C("predict", 
    as.character(Bchrondata$input[,1]),
    as.double(Bchrondata$input[,2]/1000),
    as.double(Bchrondata$input[,3]/1000),
    as.double(Bchrondata$input[,4]/100),
    as.double(Bchrondata$input[,5]/100),
    as.double(Bchrondata$input[,6]),
    as.double(Bchrondata$input[,7]),
    as.integer(Bchrondata$input[,8]),
    as.character(Bchrondata$parsfile),
    as.character(Bchrondata$chronsfile),
    as.integer(nrow(Bchrondata$input)),
    as.double(Bchrondata$outdepths/100),
    as.integer(length(Bchrondata$outdepths)),
    as.integer(numchron),
    as.double(Bchrondata$extractdate),
    as.character(Bchrondata$outlierfile))

# Finally create ranges file
# Put some of this kind of stuff in a file so that we've got depth, 2.5%, 97.5% in three columns
Bchrondata$chrons <- as.matrix(read.table(paste(Bchrondata$chronsfile)))
cat("Depth","2.5%","50%","97.5%","\n",file=paste(Bchrondata$rangesfile),append=FALSE)
for(j in 1:ncol(Bchrondata$chrons)) cat(Bchrondata$outdepths[j],quantile(Bchrondata$chrons[,j],c(0.025,0.5,0.975)),"\n",file=paste(Bchrondata$rangesfile),append=TRUE)

cat("Completed!\n")
cat("\n")

}
