\name{PlotBchron}
\alias{PlotBchron}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot chronology output from Bchron }
\description{
  Plots the chronology output from a run of \code{cpg} and \code{predict}. Best used in conjunction with the 
  \code{Bchronmenu} function for which all arguments are automatically specified.
}
\usage{
PlotBchron(FULLNAME, INFILE, DDEPTHFILE, DETSFILE, DATESFILE = NULL, RANGESFILE = NULL, COLOUR = TRUE,VERS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FULLNAME}{The full name of the core being used.}
  \item{INFILE}{ A character string detailing the location and name of the output chronologies file created 
  through \code{predict}. }
  \item{DDEPTHFILE}{ A character string detailing the location and name of the file containing the desired output
  depths. }
  \item{DETSFILE}{ A character string detailing the location and name of the ".dat" file containing the 
  radiocarbon determinations. }
  \item{DATESFILE}{ A character string detailing the location and name of the unconstrained calibrated 
  radiocarbon dates as created from \code{calibrate}. }
  \item{RANGESFILE}{ A character string detailing the location and name of the output file. }
  \item{COLOUR}{ A logical variable specifying whether the plot should be in colour or in black and white. }
  \item{VERS}{ Bchron version number used in plots. }
}

\value{
  A plot of the resultant chronology, and a file containing 4 columns; depth, 2.5th percentile, 50th percentile 
  (mode), 97.5th percentile.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }

\examples{
    ## Given in Bchronmenu    
}
\keyword{dplot}
