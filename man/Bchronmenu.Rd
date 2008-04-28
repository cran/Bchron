\name{Bchronmenu}
\alias{Bchronmenu}
\title{ Walkthrough menu system for the Bchron package }
\description{
 A walkthrough menu system for the Bchron package. Allows users the ability to load in data, run the 
 Bchron MCMC model, run predictions and plots, and look at event ages. 
}
\details{
 The Bchron menu system first creates a Bchrondata object with fields for the the input and output
 files, the name of the core being used, the details of the calibration curve, and details of the 
 model run required. Using the menu system bypasses the need to create these objects. 
} 
\usage{ Bchronmenu() } 
\author{
Andrew Parnell <Andrew.Parnell@tcd.ie>
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
