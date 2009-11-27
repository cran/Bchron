\name{Bchronmenu}
\alias{Bchronmenu}
\title{
Walkthrough menu system for the Bchron package
}
\description{
A walkthrough menu system for the Bchron package. Gives users the ability to load in data, run the 
 Bchron MCMC model, run predictions and plots, and look at event ages.
}
\usage{
Bchronmenu()
}
\details{
 The Bchron menu system first creates a Bchrondata object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. Using the menu system bypasses the need to create these objects. This function is essentially a shell for calling \code{\link{Bchronload}}, \code{\link{Bchroncalibrate}}, 
 \code{\link{BchronMCMC}}, \code{\link{Bchronpredict}} and \code{\link{Bchronplot}}. Other useful functions include: \code{\link{Bchroncheck}}, \code{\link{Bchronconvergecheck}}, 
 \code{\link{Bchroncreatecalcurve}}, \code{\link{Bchronplotevent}}, \code{\link{Bchronpredictevent}} and \code{\link{Bchronproxyplot}}
 }
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
See \code{\link{Bchronload}}, \code{\link{Bchroncalibrate}}, \code{\link{BchronMCMC}}, \code{\link{Bchronpredict}}, \code{\link{Bchronplot}}, \code{\link{Bchronpredictevent}}, 
\code{\link{Bchronplotevent}}.
}
\examples{
# For the example data, the following should run (when the files are placed in the default directory in windows; c:/Bchron)
# myrun <- Bchronload("Glendalough")
# Bchroncalibrate(myrun)
# BchronMCMC(myrun)
# Bchronpredict(myrun)
# Bchronplot(myrun)
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}

