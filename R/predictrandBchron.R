`predictrandBchron` <-
function (PARSFILE,DETSFILE,OUTFILE,lowddepths,highddepths,nddepthints,ndet,numchron,extract,OUTLIERFILE) {

cat("Parameters file at", PARSFILE, " \n")
cat("Input file at", DETSFILE, " \n")
cat("Output file at", OUTFILE, " \n")
cat("Range of design depths is \n")
for(i in 1:nddepthints) {
    cat(lowddepths[i],"to",highddepths[i],"\n")
}
cat("Outlier output file at", OUTLIERFILE, " \n")
cat("Number of determinations is", ndet, " \n")
cat("Number of design depths intervals is", nddepthints, " \n")
cat("Number of desired chronologies is", numchron, " \n")
cat("Date of extraction of core (in k yrs BP) is", extract, " \n")
    
cat("\n")

out <- .C("predictrand", 
    as.character(PARSFILE), 
    as.character(DETSFILE), 
    as.character(OUTFILE), 
    as.double(lowddepths),
    as.double(highddepths), 
    as.integer(nddepthints), 
    as.integer(ndet), 
    as.integer(numchron), 
    as.double(extract), 
    as.character(OUTLIERFILE))
}

