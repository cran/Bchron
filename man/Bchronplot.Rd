\name{Bchronplot}
\alias{Bchronplot}
\title{
Produce Bchron chronology plots
}
\description{
Produces plots of chronologies and radiocarbon dates
}
\usage{
Bchronplot(Bchrondata, plot = TRUE, dates = TRUE, datelabels = FALSE, chrons = TRUE, dateselect = NULL, limits = NULL, dateheight = mean(diff(Bchrondata$input[, 4]))/200, chronwidth = 0.95, datecolour = "black", chroncolour = "blue", outcolour = "red", transp = 0.5, outprob = 0.4, legendloc = "topleft", leghoriz = TRUE, ...)
}
\arguments{
  \item{Bchrondata}{
A Bchron data object with fields for the the input and output files, the name of the core being used, the details of the calibration curve, and details of the model run required. 
See \code{\link{Bchronload}} and \code{\link{Bchronmenu}} for more details. 
}
  \item{plot}{
If true, will create a new plot, if false, will overlay on top of existing plot
}
  \item{dates}{
If true, will plot the individually calibrated dates.
}
  \item{datelabels}{
If true, will plot the date labels to the left of the individually calibrated dates.
}
  \item{chrons}{
If true, will plot the chronology.
}
  \item{dateselect}{
A vector or sequence giving the dates which are to be plotted. Useful for ignoring the dates of boundaries, or when many dates overlap.
}
  \item{limits}{
A vector of length 4 with dimensions as c(lowestdepth,oldestage,highestdepth,youngestage).
}
  \item{dateheight}{
The height of the date densities relative to the depth scale.
}
  \item{chronwidth}{
The desired probability range of the plotted chronology. Default is 0.95 of 95%.
}
  \item{datecolour}{
The colour of the radiocarbon dates as plotted.
}
  \item{chroncolour}{
The colour of the chronology.
}
  \item{outcolour}{
The colour of outlying dates.
}
  \item{transp}{
The degree of transparency of the chronology. 
}
  \item{outprob}{
The threshold over which dates are considered outliers and coloured accordingly. Default is 0.4.
}
  \item{legendloc}{
The location of the legend. Can be given as any of "topleft", "top", etc, or locator(1) to position the legend manually.
}
  \item{leghoriz}{
If true, creates a horizontal rather than vertical legend.
}
  \item{\dots}{
Other arguments, such as width, height, of the plot created.
}
}
\details{
The options allow for a wide variety of possible plotting scenarios, including black and white plots, overlayed densities, outlier idenitfication, etc, etc.
}
\references{
Haslett and Parnell (2008, JRSSC, 57:399-418) and Parnell et al (2008, QSR, 27(19-20):1872-1885)
}
\author{
Andrew Parnell <Andrew.Parnell@ucd.ie>
}
\seealso{
See \code{\link{Bchronmenu}} for examples.
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
