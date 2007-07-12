\name{PlotDens}
\alias{PlotDens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot age density at various depths }
\description{
  Computes plots of the density for ages at various (specified) depths. Best used via \code{Bchronmenu()} which
  fills in all arguments automatically.
}
\usage{
PlotDens(datesfile, depthfile, fullname,outfile,vers,ngroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datesfile}{ A character string detailing the location and name of the output dates }
  \item{depthfile}{ A character string detailing the location and name of the desired depths }
  \item{fullname}{ The full name of the core }
  \item{outfile}{ A character string detailing the desired location and name of the outputfile }
  \item{vers}{ Bchron version number used in plots}
  \item{ngroup}{ 0 if there is only one date and you are using fixed depths, or the group number
  if you are using depth intervals and want to use a specific group}
}

\value{
  A density plot with highest posterior density ranges (HDR), as well as a file (specified in 'outfile') which
  contains HDR regions.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }
\seealso{
\code{\link{Bchronmenu}}
\code{\link{hdr2.den}}
}
\examples{
    ## Given in Bchronmenu    
}
\keyword{dplot}
\keyword{programming}
