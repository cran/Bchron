\name{BchronRSL}
\alias{BchronRSL}
\title{
Function to estimate rates of relative sea level change
}
\description{
This function takes output from a Bchron model run and fit an errors-in-variables regression of the defined order. The results can be used to estimate rates of change.
}
\usage{
BchronRSL(Bchrondata, RSLmean, RSLsd, degree = 1, iter = 10000, burnin = 2000, thin = 8, reportevery = 1000, BP = FALSE)
}
\arguments{
  \item{Bchrondata}{
The output from a Bchron model run. See \code{\link{BchronMCMC}}
}
  \item{RSLmean}{
A vector containing the mean RSL estimates. The vector should be of the same length as the number of 'output depths' of the Bchron run.
}
  \item{RSLsd}{
A vector containing the standard deviation RSL estimates. The vector should be of the same length as the number of 'output depths' of the Bchron run.

}
  \item{degree}{
The polynomial degree for the linear regression. 1 = linear, 2 = quadratic, 3 = cubic, etc
}
  \item{iter}{
The number of MCMC iterations to compute
}
  \item{burnin}{
The number of initial MCMC iterations to disgard
}
  \item{thin}{
Keep every (thin) MCMC iterations
}
  \item{reportevery}{
How often to report the results to the screen
}
  \item{BP}{
Whether to report the results in years BP (default) or years BC/AD
}
}
\details{
This function runs a non-Gaussian errors-in-variables regression on the relative sea level data taking into account the uncertainty in the chronology. The output includes error intervals for the rates, as well as the full posterior distribution in case further analysis is required.
}
\value{
A list with the following elements
  \item{samples }{The posterior samples of the regression coefficients}
  \item{degree }{The degree of the regression (see above)}
  \item{BP }{Whether run in years before present or BC/AD}
  \item{RSLmean }{The relative sea level means}
  \item{RSLsds }{The relative sea level standard deviations}
  \item{chrons }{The posterior sample of chronologies from Bchron}
  \item{const }{The mean of the posterior samples of the chronology}
}
\references{
Forthcoming
}
\author{
Andrew Parnell \email{andrew.parnell@ucd.ie}
}
\seealso{
\code{\link{BchronRSLplot}} to plot output, \code{\link{BchronMCMC}} to create Bchron chronologies.
}
\examples{
## See http://mathsci.ucd.ie/~parnell_a/
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
