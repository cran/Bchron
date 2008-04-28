\name{Bchroncalibrate}
\alias{Bchroncalibrate}
\title{ Unrestricted radiocarbon calibration of Bchron data }
\description{
  Runs standard radiocarbon calibration of determinations with no restrictions on their ordering.
}
\usage{
Bchroncalibrate(Bchrondata, iterations = 5e+05, burnin = 50000, thinby = 45, howmany = 50000)
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. }
  \item{iterations}{ Number of MCMC iterations required. }
  \item{burnin}{ Number of discarded burn-in iterations. }
  \item{thinby}{ Number of iterations to thin out. }
  \item{howmany}{ Number of iterations before R reports back and updates the screen. }
}
\references{
Haslett and Parnell (2008), A simple monotone process with application to radiocarbon-dated depth chronologies.
Applied Statistics 57 (5), 1--20.
}
\seealso{ \code{\link{BchronMCMC}} }
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
