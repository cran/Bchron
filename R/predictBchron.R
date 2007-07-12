`predictBchron` <-
function(PARSFILE,DETSFILE,OUTFILE,ndet,DDEPTHFILE,nddepths,numchron,extract,OUTLIERFILE) {

cat("Parameters file at",PARSFILE," \n")
cat("Input file at",DETSFILE," \n")
cat("Output file at",OUTFILE," \n")
cat("Pollen depths file at",DDEPTHFILE," \n")
cat("Outlier output file at",OUTLIERFILE," \n")
cat("Number of determinations is",ndet," \n")
cat("Number of pollen depths is",nddepths," \n")
cat("Number of desired chronologies is",numchron," \n")
cat("Date of extraction of core (in k yrs BP) is",extract," \n")

cat("\n")

out <- .C("predict",
    as.character(PARSFILE),
    as.character(DETSFILE),
    as.character(OUTFILE),
    as.integer(ndet),
    as.character(DDEPTHFILE),
    as.integer(nddepths),
    as.integer(numchron),
    as.double(extract),
    as.character(OUTLIERFILE))
}

