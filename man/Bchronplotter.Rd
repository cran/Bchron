\name{Bchronplotter}
\alias{Bchronplotter}
\title{ Produce Bchron chronology plots }
\description{
Produces colour and B/W plots of a chronology
}
\usage{
Bchronplotter(Bchrondata, colours = c("yellow", "orange", "red"), lwidths = c(3, 3, 3), asklimits = FALSE)
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. }
  \item{colours}{ Colours used for the 95, 75 and 50 percent highest posterior density regions (HDRs). }
  \item{lwidths}{ Widths of the lines used for the 95, 75 and 50 percent HDRs. }
  \item{asklimits}{ Whether to request plotting ranges in terms of depths and ages. }
}
\references{
Haslett and Parnell (2008), A simple monotone process with application to radiocarbon-dated depth chronologies.
Applied Statistics 57 (5), 1--20.
}
\seealso{ \code{\link{Bchronpredict}} }
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
