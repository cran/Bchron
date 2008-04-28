\name{Bchronconvergecheck}
\alias{Bchronconvergecheck}
\title{ Convergence checking of Bchron output }
\description{
  Allows convergence checking of both \code{\link{BchronMCMC}} and \code{\link{Bchroncalibrate}} output. Uses the
  Geweke convergence diagnostic from the coda package.
}
\usage{
Bchronconvergecheck(Bchrondata = list(RUN = FALSE), CALDATES = FALSE)
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required.  }
  \item{CALDATES}{ TRUE if the calibrated dates are being checked (after a run of \code{\link{Bchroncalibrate}}), 
  FALSE if a check of \code{\link{BchronMCMC}} is required. }
}
\references{
Haslett and Parnell (2008), A simple monotone process with application to radiocarbon-dated depth chronologies.
Applied Statistics 57 (5), 1--20.
}
\seealso{ \code{\link{Bchronloaddata}} }
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}

