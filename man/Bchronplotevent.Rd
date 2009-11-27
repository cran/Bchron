\name{Bchronplotevent}
\alias{Bchronplotevent}
\title{
Histogram of the probability distribution for an event in a core.
}
\description{
Provides a histogram and summary of the probability distribution for a given depth, slice, or event in a core.
}
\usage{
Bchronplotevent(Bchrondata, depth = NULL, slice = NULL, event = NULL, eventname = NULL, nbreaks = 30, histcolour = "light blue", ...)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details. 
}
  \item{depth}{
One of the output depths (in cm) in the core (see \code{\link{Bchronload}} for more information on output depths).
}
  \item{slice}{
A slice number for one of the depths in the core.
}
  \item{event}{
An event name. Must have an associated corenameeventdepthseventname.txt file associated with it.
}
  \item{eventname}{
The full name of the event (eg "Alnus rise").
}
  \item{nbreaks}{
The number of breaks for the histogram.
}
  \item{histcolour}{
The colour of the histogram.
}
  \item{\dots}{
Other arguments, such as the width and height of the plot.
}
}
\details{
This function can be used to identify the age of a certain part of the core according to the chronology model. It does not necessarily have to be an age that has a radiocarbon (or non-radiocarbon) 
date. Where an event name is given, a file should be created giving the depths at which the event occurred. See \code{\link{Bchronpredictevent}} for more information.
}
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
\code{\link{Bchronpredictevent}}
}
\examples{
# For the example data, the following should run (when the files are placed in the default directory in windows; c:/Bchron)
# myrun <- Bchronload("Glendalough")
# Bchroncalibrate(myrun)
# BchronMCMC(myrun)
# Bchronpredictevent(myrun,event="Alnus")
# Bchronplotevent(myrun,event="Alnus",eventname="Alnus rise")
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}

