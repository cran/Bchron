\name{BchronMCMC}
\alias{BchronMCMC}
\title{
Run the Bchron
}
\description{
 Runs the compound Poisson-Gamma model for radiocarbon and non-radiocarbon dates
}
\usage{
BchronMCMC(Bchrondata, iterations = 1e+05, burnin = 10000, thinby = 8, howmany = 2000, defaults = FALSE, testrun = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details.
}
  \item{iterations}{
    Number of MCMC iterations required. Default is 100,000.
}
  \item{burnin}{
    Number of discarded burn-in iterations. Default is 10,000.
}
  \item{thinby}{
    Step size of iterations to thin out. Default is to take every 8th iteration.
}
  \item{howmany}{
    Number of iterations before R reports back and updates the screen. Default is 2000.
}
  \item{defaults}{
    If TRUE, runs BchronMCMC without any warnings and overwrites previous runs. If FALSE (default) will prompt user before continuing.
}
  \item{testrun}{
    if TRUE, tries a very short run to test that BchronMCMC works.
}
}
\details{
    The Bchron MCMC function fits a compound Poisson-Gamma distribution to the increments between the dated levels. 
    This involves a stochastic linear interpolation where the age gaps are Gamma(2,beta) distributed, and the 
    depth gaps are Exponential(lambda). Radiocarbon and non-radiocarbon dates (including outliers) are updated within the function also by MCMC.
}
\value{
Creates a parameter output file (corenamepars.txt) in the output directory containing outlier estimates, parameter values, and other details required. This information can be turned into chronology
samples by calling \code{\link{Bchronpredict}}.
}
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
See \code{\link{Bchronmenu}} for examples.
}
\examples{
# For the example data, the following should run (when the files are placed in the default directory in windows; c:/Bchron)
# myrun <- Bchronload("Glendalough")
# Bchroncalibrate(myrun)
# BchronMCMC(myrun)
# Bchronpredict(myrun)
# Bchronplot(myrun)
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
