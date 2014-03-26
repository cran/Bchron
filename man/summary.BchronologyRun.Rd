\name{summary.BchronologyRun}
\alias{summary.BchronologyRun}

\title{
Summarise a Bchronology object
}
\description{
Summarise a \code{\link{Bchronology}} object
}
\usage{
\method{summary}{BchronologyRun}(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), 
\ldots, digits = max(3, getOption("digits") - 3))
}
\arguments{
  \item{object}{
Output from a run of \code{\link{Bchronology}}
}
  \item{probs}{
(optional) Probabilities (between 0 and 1) at which to summarise the predicted chronologies
}
  \item{\dots}{
Other arguments (not currently supported)
}
  \item{digits}{
Number of digits to report values (not currently supported)
}
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie
}
\seealso{
\code{\link{BchronCalibrate}}, \code{\link{Bchronology}} \code{\link{BchronRSL}}, \code{\link{BchronDensity}}, \code{\link{BchronDensityFast}} 
}
\examples{
# Examples in main function \code{\link{Bchronology}}
}
