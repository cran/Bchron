\name{Bchronpredictevent}
\alias{Bchronpredictevent}
\title{
Predict ages of events in a core 
}
\description{
 Estimates and plots the ages of specific depths, slices or events in a Bchron core. Useful for event identification.
}
\usage{
Bchronpredictevent(Bchrondata, event = "Event", numchron = 10000)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details. 
}
  \item{event}{
An event name. Must have an associated corenameeventdepthseventname.txt file associated with it.
}
  \item{numchron}{
The number of chronologies required. Default is 10,000.
}
}
\details{
For each event, the user will need to create a file of the desired depths intervals for which the ages are required. These depths are to be put the files in the input directory in a 
file named CorenameEventDepthsEventname.txt where Corename is the name of the core and Eventname is the name of the event (eg GlendaloughEventDepthsUlmusDecline.txt). The format of the 
file should look like the following:
261 265
267 269
271 282
This formulation allows for uncertainty in the depths at which the event occurred.
}
\value{
A file CorenameEventAgesEventName.txt will be created in the output directory giving extra information on the ages of the event.
}
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
\code{\link{Bchronplotevent}}
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
