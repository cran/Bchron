print.Bchron <-
function(x,...){
  cat('Core name is',x$fullname, ' \n')
  cat("Calibration curve at", x$calibcurvefile, " \n")
  cat("Input file at", x$inputfile, " \n")
  cat("Output file at", x$parsfile, " \n")
  cat("Number of determinations is", nrow(x$input), " \n")
  cat("\n") 
}
