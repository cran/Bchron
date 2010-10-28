\name{Bchron-package}
\alias{Bchron-package}
\alias{Bchron}
\docType{package}
\title{
Compound Poisson-Gamma chronologies for radiocarbon and non-radiocarbon dates.
}
\description{
Runs MCMC and predictive distributions for the monotonic piece-wise linear compound Poisson Gamma process.
Before a run of the Bchron model can commence, the file structures must be in the correct format. Included in the package
is an example set of data from Glendalough, Ireland. All other cores must be in a similar format for the routines
to work. The latest version of Bchron also contains a function (\code{\link{Bchronproxyplot}}) for plotting proxy data with age uncertainty.

 ***IMPORTANT NOTE*** This version fixes a serious bug which was present in versions 3.1.1 and 3.1.2. You should have noticed this bug anyway when you ran the Bchronconvergecheck function.
}
\details{
\tabular{ll}{
Package: \tab Bchron\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2009-11-27\cr
License: \tab GPL (>2) \cr
LazyLoad: \tab yes\cr
}
Bchron will run on windows, mac or linux systems, though installation instructions are slightly different for each. Instructions for completing installation are as follows:

1. Create a directory on your hard drive called Bchron. 

2. Now navigate to the R directory and find the Bchron sub-directory within. On windows this will be C:\\ program files\\ R\\ R-XXXXX\\ library\\ BChron\\ where XXXXX is the 
version number of R. On other plaforms, this directory can be found by typing .libPaths() at the R command prompt.

3. In this directory, there there should be three subdirectories called Input, Output, and CalCurve. Copy these to the Bchron directory you created in step 1.

Everything is now set up.

Example full model run (with Glendalough.dat):

1. At the command prompt in R, type library(Bchron)

2. Type Bchronmenu() and choose option 1.

3. Locate your Bchron directory, and select IntCal04.bch as your calibration curve and Glendalough.dat as your input file. All other options
can be left as default.

4. Now choose option 2 to calibrate the radiocarbon dates.

5. Choose option 3 and 'short' to do a short run of the Bchron model.

6. Choose option 4 to run the prediction stage and create a plot of the data.

Once a satisfactory short run has been obtained, a long run should be undertaken. The long run will take much longer
than the short run, but will only be required once.

Example event prediction stage (with GlendaloughEventDepthsAlnus.txt)

1. Follow the steps above, and now choose option 5 for event prediction.

2. Type Alnus as the name of the event.

3. Check the output directory for GlendaloughEventAgesAlnus.txt which will contain 95\% HDR age intervals for the depths of interest.

Other cores can be created in a spirit similar to that found in the Glendalough.dat file. This example input file has
5 radiocarbon dates and a non-radiocarbon date for the top of the core. The columns are tab
delimited and represent the lab code of the sample, the radiocarbon age, the sample standard error, the depth (in cm)
at which it was found, the thickness of the sample in cm. If the thickness is unknown, zero is acceptable. The next
two columns contain the probabilities of being an outlier. The first probability identifies censored outliers as proposed
by Christen (1994), the second indicates the probability an outlier is ignored completely by the Bchron model. The final
column indicates whether this date is a radiocarbon date (type 1), a uniformly distributed date (type 2), or a normally
distributed date (type 3). For uniform dates, the standard deviation value is taken to be the distance to the upper and 
lower limits. Type 3 dates are permitted to be outliers, with probabilities as given. Type 2 dates are not allowed to be
outliers; the outlier probabilities for these dates are ignored.
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
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}


