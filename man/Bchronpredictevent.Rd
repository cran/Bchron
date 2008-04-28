\name{Bchronpredictevent}
\alias{Bchronpredictevent}
\title{ Predict the ages of events in Bchron }
\description{
    Estimates the ages of specific depths in a Bchron core. Useful for event identification.
}
\usage{
Bchronpredictevent(Bchrondata)
}
\arguments{
  \item{Bchrondata}{ A Bchron data object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. }
}
\details{
  Produces a file AeventagesB.txt where A is the name of the core and B is the name of the event. File 
  contains sampled ages for this event. Note that the depths of the event in the core can be specified in 
  terms of ranges for scenarios where the event is not easily identified. 
}
\references{
Haslett and Parnell (2008), A simple monotone process with application to radiocarbon-dated depth chronologies.
Applied Statistics 57 (5), 1--20.
}
\seealso{ \code{\link{Bchronplotevent}} }
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
