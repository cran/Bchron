\name{predictBchron}
\alias{predictBchron}
\alias{predict}
\title{ Predict calendar ages at desired fixed depths }
\description{
Best accessed through the \code{Bchronmenu} where arguments are filled in automatically. This function will take a set of parameter values created via \code{Bchron}
and draw stochastically consistent chronologies. Note that this function is generally pretty fast, runs will usually
take less than a minute. 
}
\usage{
predictBchron(PARSFILE, DETSFILE, OUTFILE, ndet, DDEPTHFILE, nddepths, numchron, extract,OUTLIERFILE)
}
\arguments{
  \item{PARSFILE}{ A character string detailing the location and name of the parameters file from a run of the Bchron model. }
  \item{DETSFILE}{ A character string detailing the location and name of the determinations file consistent with that 
  expected by the Bchron model. See Bchronmenu for more details. }
  \item{OUTFILE}{ A character string detailing the desired location and name of the output chronologies. }
  \item{ndet}{ The number of determinations }
  \item{DDEPTHFILE}{ A character string detailing the location and name of the file containing the desired depths at which
  to draw the chronology. }
  \item{nddepths}{ The number of desired depths }
  \item{numchron}{ The number of chronologies desired. }
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
