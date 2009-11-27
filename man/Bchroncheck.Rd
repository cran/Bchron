\name{Bchroncheck}
\alias{Bchroncheck}
\title{
Checking of Bchron data objects
}
\description{
    A quick check to make sure that the Bchron data object will run in \code{\link{Bchroncalibrate}} and \code{\link{BchronMCMC}}. 
}
\usage{
Bchroncheck(Bchrondata)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details.
}
}
\details{
Checks the Bchron data object for inconsistencies. If the object passes all these tests, it *should* run correctly in all areas of the Bchron system.
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

