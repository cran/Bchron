\name{Bchronproxyplot}
\alias{Bchronproxyplot}
\title{
Plot proxy data with associated chronological unceratinties
}
\description{
Function produces plots of a proxy data time series, with associated error bands to indicate chronological uncertainty.
}
\usage{
Bchronproxyplot(Bchrondata, proxy, title = NULL, xlabel = "Age (k cal yrs BP)", ylabel = "Proxy", num = 10, col = c("red","red", "red"), smooth = FALSE,...)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details. 
}
  \item{proxy}{
A vector of proxy observations of the same length as the number of output depths in the core. See \code{\link{Bchronload}} for more information on the output depths.
}
  \item{title}{
A text string giving the title for the plot.
}
  \item{xlabel}{
A text string giving the x-axis label for the plot
}
  \item{ylabel}{
A text string giving the y-axis label for the plot
}
  \item{num}{
The number of sample chronologies to plot
}
  \item{col}{
The colour of the lines drawn
}
  \item{smooth}{
Whether kernel smoothing on the proxy variable should be performed
}
  \item{\dots}{
Other arguments, such as the width and height of the plot
}
}
\details{
The plots produces show a number of sample proxy series based on different chronological samples (the number of these is controlled by num), together with an estimate of the 95% range of 
uncertainty in the proxy series. As the uncertainty is in the horizontal axis, uncertainty bounds can sometimes cross each other if there is an inflexion point.
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
# x <- sin(seq(0,pi,length=150))
# Bchronproxyplot(myrun,x)
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}

