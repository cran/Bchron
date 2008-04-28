\name{Bchronplotterdatesonly}
\alias{Bchronplotterdatesonly}
\title{ Produce Bchron plots of calibrated radiocarbon dates}
\description{
Produces colour and B/W plots of the independently calibrated radiocarbon dates.
}
\usage{
Bchronplotterdatesonly(Bchrondata, asklimits = FALSE)
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. }
  \item{asklimits}{ Whether to request plotting ranges in terms of depths and ages. }
}
\references{
Haslett and Parnell (2008), A simple monotone process with application to radiocarbon-dated depth chronologies.
Applied Statistics 57 (5), 1--20.
}
\seealso{ \code{\link{Bchronplotter}} }
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
