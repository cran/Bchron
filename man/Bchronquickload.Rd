\name{Bchronquickload}
\alias{Bchronquickload}
\title{ Quick loading of Bchron data }
\description{
    A shortcut to \code{\link{Bchronloaddata}} for advanced users.
}
\usage{
Bchronquickload(name, fullname = NULL, path = "c:/Bchron", calibname = "IntCal04", Bchronversion = 0, check = FALSE)
}
\arguments{
  \item{name}{ Name of the files which contain the input data details (eg "Glendalough") }
  \item{fullname}{ Full name fo the core }
  \item{path}{ The path to the input, output and calcurve directories. }
  \item{calibname}{ The name of the calibration curve. }
  \item{Bchronversion}{ The version number of Bchron (optional; used for plotting only). }
  \item{check}{ Whether data checks should be performed on the data. }
}
\details{
    Output should be analagous to the \code{\link{Bchronloaddata}} function. 
}
\value{
  Returns a Bchrondata object with some of the following
  \item{SHOULDRUN}{(logical) Determines whether the \code{\link{BchronMCMC}} can be run.}
  \item{path }{(character) The path to the input, output and calcurve directories. }
  \item{calibcurvefile }{(character) The path to the calibration curve file. }
  \item{bigcalsize }{The number of entries in the calibration curve file. }
  \item{lowcal }{The lowest calendar age limit allowed on the calibration curve in k cal years BP (-0.005 for IntCal04). }
  \item{highcal }{The highest calendar age limit allowed on the calibration curve in k cal years BP (26 for IntCal04). }
  \item{name }{(character) The name of the input files for the core. }
  \item{ddepthfile }{(character) The path to the design depths file (the dates at which pollen was taken for a pollen core). }
  \item{inputfile }{(character) The path to the input file where the 14C dates are stored. }
  \item{fullname }{(character) The full name of the core. }
  \item{calibdatesfile }{(character) The path to the unrestricted calibrated dates file. }
  \item{parsfile }{(character) The path to the output file where the parameters are/will be stored. }
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
