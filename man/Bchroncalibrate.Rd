\name{Bchroncalibrate}
\alias{Bchroncalibrate}
\title{
Radiocarbon calibration of Bchron data
}
\description{
Runs standard radiocarbon calibration of determinations.
}
\usage{
Bchroncalibrate(Bchrondata, iterations = 5e+05, burnin = 50000, thinby = 45, howmany = 50000, defaults = FALSE)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details. 
}
  \item{iterations}{
    Number of MCMC iterations required. Default is 500,000.
}
  \item{burnin}{
    Number of discarded burn-in iterations. Default is 50,000.
}
  \item{thinby}{
    Step size of iterations to thin out. Default is to take every 45th iteration.
}
  \item{howmany}{
    Number of iterations before R reports back and updates the screen. Default is 50000
}
  \item{defaults}{
    If true, runs Bchroncalibrate without any warnings and overwrites previous runs. If FALSE (default) will prompt user before continuing.
}
}
\details{
This function run simple radiocarbon calibration of each of the determinations individually without any chronology constraint. The results here should be directly comparable with those 
obtained from, eg Oxcal, Bcal, Calib, etc.
}
\value{
No specific value is returned but two files are created in the output directory, giving the raw radiocarbon calibrated age iterations (corenameTrueDates.txt) and a summary of the 
2.5%, 50%, and 97.5% probability ranges for each date (corenameCalibRanges.txt).
}
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for examples.
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
