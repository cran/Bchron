\name{Bchronchangecalcurve}
\alias{Bchronchangecalcurve}
\title{ Change Bchron radiocarbon calibration curve }
\description{
    Allows the user to switch between calibration curves.
}
\usage{
Bchronchangecalcurve(Bchrondata)
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. }
}
\details{
  The default calibration curve is IntCal04. Bchron uses a specific, yearly-interpolated version of this calibration
  curve to speed up computation. This function allows users to download other internationally-agreed calibration 
  curves (as .14c files) and convert them into Bchron .bch format. The new .14c files must be put in the Bchron 
  calcurve directory. 
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
