\name{hdr2.den}
\alias{hdr2.den}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Highest density region plots for Bchron output. }
\description{
  Adapted from \code{hdr.den} in the \code{hdrcde} package which is a pre-requisite for the Bchron package
  to work. Produces plots of the highest density regions for ages at desired depths throughout the core.
}
\usage{
hdr2.den(x = NULL, prob = c(50, 95, 99), den = NULL, h = NULL, ylab = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Simulations from the density of interest. }
  \item{prob}{ Pobability regions }
  \item{den}{ The density for \code{x}. If not given, is automatically calculated }
  \item{h}{ Bandwidth }
  \item{ylab}{ The y-axis label }
  \item{\dots}{ Other arguments passed to the \code{plot} command. }
}
\value{
  A plot of the density for together with 50\%, 95\% and 99\% highest posterior density regions.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }

\examples{
    ## Given in Bchronmenu    
}
\keyword{dplot}
\keyword{nonparametric}
\keyword{smooth}
