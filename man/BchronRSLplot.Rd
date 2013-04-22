\name{BchronRSLplot}
\alias{BchronRSLplot}
\title{
Plot of Bchron RSL output
}
\description{
Plots of the output from a run of \code{\link{BchronRSL}}. 
}
\usage{
BchronRSLplot(Bchrondata, RSLresults, xgrid = NULL, colours = colorRampPalette(c("white", "black", "black", "black", "black")), transp = 0.2, ss = 1000, ellipses = TRUE, blobs = FALSE)
}
\arguments{
  \item{Bchrondata}{
The output from a Bchron model run. See \code{\link{BchronMCMC}}
}
  \item{RSLresults}{
The output from a BchronSL model run. See \code{\link{BchronRSL}}
}
  \item{xgrid}{
A grid of ages on which to interpolate and produce predictions
}
  \item{colours}{
The colours used to draw the blobs
}
  \item{transp}{
The degree of transparancy of the blobs
}
  \item{ss}{
The sample size used to create the kernel density esimates of the blobs
}
  \item{ellipses}{
Whether to draw uncertainty ellipses around each RSL point
}
  \item{blobs}{
Whether to draw the full probability distribution for each RSL point. This can be slow if set to TRUE
}
}
\value{
No output produced, just a plot
}
\references{
Forthcoming
}
\author{
Andrew Parnell \email{andrew.parnell@ucd.ie}
}
\seealso{
\code{\link{BchronRSLplot}} to plot output, \code{\link{BchronMCMC}} to create Bchron chronologies.
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
