\name{BchronreadEPD}
\alias{BchronreadEPD}
\title{ Reading in of European pollen database (EPD) data }
\description{
  A menu driven function that creates Bchron-style input files from European Pollen Database (EPD) data. 
  Note that the EPD data must be saved to a file containing only the dating information.
}
\usage{
BchronreadEPD()
}
\arguments{
}
\details{
  This function allows the user to cut or save data direct from the EPD and then transform it into a Bchron-style input 
  file. Note that some of the features of the data (such as top dates, or non-radiocarbon dates) may not be preserved, 
  so manual editing of these files may be necessary.
}
\value{
    Creates a .dat file (as specified by the user) which contains the radiocarbon dates and depths.
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
