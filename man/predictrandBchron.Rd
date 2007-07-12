\name{predictrandBchron}
\alias{predictrandBchron}
\title{ Predict calendar ages at desired random uniform interval depths }
\description{
Best accessed through the \code{Bchronmenu} where arguments are filled in automatically. This function will take a 
set of parameter values created via \code{Bchron} and draw stochastically consistent chronologies. Note that this 
function is generally pretty fast, runs will usually take less than a minute. In particular, this version takes
low and high depths (in a list, there can be more than one) and generates design depths which are uniformly 
distributed between these values. This is used where an event can only be identified as occurring between two
depths but not at an exact one.
}
\usage{
predictrandBchron(PARSFILE, DETSFILE, OUTFILE, lowddepths, highddepths, nddepthints, ndet, numchron, extract, OUTLIERFILE)
}
\arguments{
  \item{PARSFILE}{ A character string detailing the location and name of the parameters file from a run of the Bchron model. }
  \item{DETSFILE}{ A character string detailing the location and name of the determinations file consistent with that 
  expected by the Bchron model. See Bchronmenu for more details. }
  \item{OUTFILE}{ A character string detailing the desired location and name of the output chronologies. }
  \item{lowddepths}{ A list containing the lower bounds for the event depth. }
  \item{highddepths}{ A list containing the upper bounds for the event depth. }
  \item{nddepthints}{ The number of different intervals for the event depth. }
  \item{ndet}{ The number of determinations }
  \item{numchron}{ The number of chronologies desired.  }
  \item{extract}{ The date of extraction of the core. This is used to stop chronologies being drawn for impossible ages 
  in the future. The procedure may be slowed down somewhat if the model deems that these types of ages are probable. }
  \item{OUTLIERFILE}{ A character string detailing the desired location and name of the outlier summary of the core. }
}
\value{
  No output produced by the function, but will populate the file designated OUTFILE and OUTLIERFILE.
}
\references{ Haslett and Parnell (2007, submitted to JRSSC) }
\author{ Andrew Parnell }

\seealso{ \code{\link{Bchronmenu}} }
\examples{    
    ## Given in Bchronmenu
}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
