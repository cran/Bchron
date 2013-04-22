\name{Quickcal}
\alias{Quickcal}
\title{
Quick radiocarbon calibration function
}
\description{
Runs a quick radiocarbon calibration estimate for a single date and standrad deviation.
}
\usage{
Quickcal(date, sd, pathtocalcurve=NULL,prob=c(68,95,99),get.output=FALSE)
}
\arguments{
  \item{date}{
The 14C age of the sample
}
  \item{sd}{
The 1-sigma standard deviation (in 14C years) of a sample
}
  \item{pathtocalcurve}{
A character string giving the file path to the radiocarbon calibration curve
}
  \item{prob}{
The values of the highester posterior density regions to compute for the calibrated date
}
  \item{get.output}{
Whether to provide the date density as an output
}
}
\details{
This function provides a very fast radiocarbon calibration for individual dates and standard deviations given in radiocarbon years for a given calibration curve. Suitable .14c calibration curves can be found on the web, though the latest IntCal northern hesmispher curve is included as part of the Bchron package. See \code{\link{BchronMCMC}} for more information on how to obtain this calibration curve.
}
\value{
If get.output=TRUE, then a list with the following elements
  \item{age }{THe range of ages with positive probability density}
  \item{density }{The value of the probability density for the given age}
}
\references{
Forthcoming
}
\author{
Andrew Parnell \email{andrew.parnell@ucd.ie}
}
\seealso{
\code{\link{BchronMCMC}} to for details on accessing the calibration curve.
}
\examples{
## See http://mathsci.ucd.ie/~parnell_a/
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
