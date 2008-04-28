\name{Bchronpredict}
\alias{Bchronpredict}
\title{ Predict ages for depths in a Bchron core }
\description{
    Performs stochastic Poisson-Gamma linear interpolation on the output of a run of the \code{\link{BchronMCMC}}.
}
\usage{
Bchronpredict(Bchrondata = list(SHOULDRUN = FALSE))
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. }
}
\details{
    Produces a file *chrons.txt which contains sample chronologies for each of the design depths.
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
