\name{Bchronconvergecheck}
\alias{Bchronconvergecheck}
\title{
MCMC convergence checking of Bchron output
}
\description{
Allows convergence checking of both \code{\link{BchronMCMC}} and \code{\link{Bchroncalibrate}} output. Uses the
Geweke convergence diagnostic from the coda package.
}
\usage{
Bchronconvergecheck(Bchrondata, dates = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details.
}
  \item{dates}{
if true, checks the individually calibrated dates for convergence. Otherwise checks the parameter set from a run of \code{\link{BchronMCMC}}.
}
}
\details{
Uses the package \code{coda} Geweke diagnostic to perform a convergence check on a single MCMC chain. Note that this is not the only method of convergence checking. Also provides some
context for how well the model appears to have converged.
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


