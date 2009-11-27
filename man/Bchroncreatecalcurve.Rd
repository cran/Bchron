\name{Bchroncreatecalcurve}
\alias{Bchroncreatecalcurve}
\title{
Create Bchron-style radiocarbon calibration curve}
\description{
Allows the user to create calibration curves for use with Bchron.
}
\usage{
Bchroncreatecalcurve()
}
\details{
This function will turn a standard .14c file (as given by the international radiocarbon calibration group) into a .bch file as used by Bchron. Bchron uses a specific, 
yearly-interpolated version of the calibration curves to speed up computation. Bchron comes with northern-hemisphere, southern-hemisphere, and marine reservoir calibration curves, but users 
can create their own calibration files in the .14c format if they wish. This function will translate that calibration curve into the Bchron .bch format.
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
