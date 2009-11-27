\name{Bchronpredict}
\alias{Bchronpredict}
\title{
Predict ages for given depths in a Bchron core
}
\description{
Performs stochastic Poisson-Gamma linear interpolation on the output of a run of the \code{\link{BchronMCMC}}.
}
\usage{
Bchronpredict(Bchrondata, numchron = 10000, defaults = FALSE)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details. 
}
  \item{numchron}{
The number of chronologies required. Default is 10,000.
}
  \item{defaults}{
    If true, runs Bchronpredict without any warnings and overwrites previous runs. If FALSE (default) will prompt user before continuing.
}
}
\details{
Must be run after a \code{\link{BchronMCMC}} run to turn parameter estimates into sample chronologies.    
}
\value{
Produces files corenamechrons.txt and corenameranges.txt in the output directory which contain sample chronologies and estimated ranges for each of the output depths.
}
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
\code{\link{Bchronpredictevent}}
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
