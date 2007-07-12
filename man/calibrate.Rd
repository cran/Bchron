\name{calibrate}
\alias{calibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calibrate a set of radiocarbon determinations }
\description{
  The calibrate function will take a set of radiocarbon determinations and calibrate them via MCMC. The output
  should be *consistent* with that produced by other calibration methods, eg bcal, oxcal, calib. 
}
\usage{
calibrate(CALPATH, INFILE, OUTFILE, ndet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CALPATH}{ A character string detailing the location and name of the calibration curve file.
  Note that Bchron uses a specific version of the Northern hemisphere Intcal calibration curve. Other calibration
  curves are not supported. }
  \item{INFILE}{ A character string detailing the location and name of the determinations, as required by
  Bchron package. See \code{\link{Bchronmenu}} for more details. }
  \item{OUTFILE}{ A character string detailing the desired location and name of the output file. }
  \item{ndet}{ The number of radiocarbon determinations for this core. }
}

\value{
   No output produced but the file specified in \code{OUTFILE} is populated by calibrated dates.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }

\examples{
    ## Given in Bchronmenu    
}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
