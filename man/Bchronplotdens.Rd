\name{Bchronplotdens}
\alias{Bchronplotdens}
\title{ Event density plotting for Bchron data }
\description{
    Not meant for use other than via \code{\link{Bchronplotevent}}
}
\usage{
Bchronplotdens(datesfile, fullname, outfile, eventname, vers = 0)
}
\arguments{
  \item{datesfile}{ Character string containing file and path to ages of event. }
  \item{fullname}{ The full name of the core (used for plotting). }
  \item{outfile}{ Character string specifying path and name of file for highest posterior density regions of 
    the age of the event. }
  \item{eventname}{ Name of the event. }
  \item{vers}{ Bchron version number (used for plotting). }
}
\references{
Haslett and Parnell (2008), A simple monotone process with application to radiocarbon-dated depth chronologies.
Applied Statistics 57 (5), 1--20.
}
\seealso{ \code{\link{Bchronloaddata}}, \code{\link{Bchronquickload}}}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
