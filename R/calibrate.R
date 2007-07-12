`calibrate` <-
function(CALPATH,INFILE,OUTFILE,ndet) {

cat("Calibration curve at",CALPATH," \n")
cat("Input file at",INFILE," \n")
cat("Output file at",OUTFILE," \n")
cat("Number of determinations is",ndet," \n")
cat("\n")

out <- .C("calibrate",
    as.character(CALPATH),
    as.character(INFILE),
    as.character(OUTFILE),
    as.integer(ndet))
}

