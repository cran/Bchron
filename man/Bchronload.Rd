\name{Bchronload}
\alias{Bchronload}
\title{
Command-line based loading of Bchron input files for a Bchron run.
}
\description{
A quick command line function to create a Bchron data object from various inputs. Note that all text strings must be given in inverted commas (apart from TRUE/FALSE).
}
\usage{
Bchronload(name, fullname = NULL, path = NULL, outdepths=NULL,calibname = "IntCal09", extractdate = -0.05, check = FALSE, full = FALSE)
}
\arguments{
  \item{name}{
A text string containing the name of the .dat file used.
}
  \item{fullname}{
A text string giving the full name of the core (used for plotting purposes). Defaults to name.
}
  \item{path}{
The path to the Bchron directory. Defaults to c:/Bchron (windows) or ~/Bchron (unix). 
}
  \item{outdepths}{
A collection of depths at which chronological ages will be calculated. Default is top dated depth to bottom dated depth, length 200.
}
  \item{calibname}{
The name of the calibratiion curve .bch file used. Defaults to the Northern Hemisphere 2009 internationally agreed calibration curve.
}
  \item{extractdate}{
The extracted date of the core in thousands of years before present. Defaults to the year 2000 (=-0.05 k cal years BP).
}
  \item{check}{
If true, performs \code{\link{Bchronconvergecheck}} on the given data files.
}
  \item{full}{
If true, also loads in chronologies, calibrated dates, etc, for further analysis. This may take a bit of extra time as some of these files are large.
}

}
\details{
Useful for batch processing of calibration runs, this function allows the user to load in and create a Bchrondata object which can then be used in all the other Bchron functions.
}
\value{
Returns a list object with the following components:
  \item{name}{The name of the .dat file}
  \item{fullname}{The full name of the core}
  \item{version}{The version number of Bchron}
  \item{path}{The path to the Bchron directory}
  \item{calibcurvefile}{The path to the calibration curve}
  \item{c14file}{The path to the .14c file of the calibration curve}
  \item{inputfile}{The path to the input file of dates and depths}
  \item{calibdatesfile}{The path to the file of calibrated radiocarbon date samples}
  \item{calibrangesfile}{The path to the file of calibrated radiocarbon date ranges}
  \item{parsfile}{The path to the file of Bchron parameter output}
  \item{chronsfile}{The path to the file containing the sample chronologies}
  \item{outlierfile}{The path to the file containing the outlier date sumaries}
  \item{rangesfile}{The path to the file containing the chronology summary ranges}
  \item{extractdate}{The date of core extraction}
  \item{bigcalsize}{The length of the .bch calibration curve}
  \item{lowcal}{The youngest age of the given calibration curve}
  \item{highcal}{The eldest age of the given calibration curve}
  \item{input}{The table of input data (dates, depths, thicknesses etc)}
  \item{outdepths}{A collection of depths at which the chronological age is required}
  \item{calibdates}{(if full=TRUE, a matrix of sample calibrated dates}
  \item{calibranges}{if full=TRUE, a matrix of calibrated date ranges}
  \item{pars}{if full=TRUE, a matrix of Bchron parameter output}
  \item{chrons}{if full=TRUE, a matrix of sample chronologies}
  \item{outlier}{if full=TRUE, a matrix of outlier date summaries}
  \item{ranges}{if full=TRUE, a matrix of chronology summaries}
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
