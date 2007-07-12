`BchronMCMC` <-
function(CALPATH,INFILE,OUTFILE,ndet,iterations,burnin,howmany,thinby) {

cat("Calibration curve at",CALPATH," \n")
cat("Input file at",INFILE," \n")
cat("Output file at",OUTFILE," \n")
cat("Number of determinations is",ndet," \n")
cat("\n")

out <- .C("cpg",
    as.character(CALPATH),
    as.character(INFILE),
    as.character(OUTFILE),
    as.integer(ndet),
    as.integer(iterations),
    as.integer(burnin),
    as.integer(howmany),
    as.integer(thinby))
}

