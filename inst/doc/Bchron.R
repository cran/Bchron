## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height = 5, fig.align = 'center',
                      warning=FALSE, message=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  install.packages('Bchron')

## ------------------------------------------------------------------------
library(Bchron)

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
ages1 = BchronCalibrate(ages=11553,
                        ageSds=230,
                        calCurves='intcal13',
                        ids='Date-1')
summary(ages1)
plot(ages1)

## ----results='hide'------------------------------------------------------
ages2 = BchronCalibrate(ages=c(3445,11553,7456), 
                        ageSds=c(50,230,110), 
                        calCurves=c('intcal13','intcal13','shcal13'))
summary(ages2)
plot(ages2)

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
ages3 = BchronCalibrate(ages=c(3445,11553), 
                        ageSds=c(50,230), 
                        positions=c(100,150), 
                        calCurves=c('intcal13','normal'))
summary(ages3)
plot(ages3,withPositions=TRUE)

## ------------------------------------------------------------------------
# First create age samples for each date
age_samples = SampleAges(ages3)
# Now summarise them with quantile - this gives a 95% credible interval
apply(age_samples,2,quantile,prob=c(0.025,0.975))

## ------------------------------------------------------------------------
data(Glendalough)
print(Glendalough)

## ----results='hide'------------------------------------------------------
GlenOut = Bchronology(ages=Glendalough$ages,
                      ageSds=Glendalough$ageSds, 
                      calCurves=Glendalough$calCurves,
                      positions=Glendalough$position, 
                      positionThicknesses=Glendalough$thickness,
                      ids=Glendalough$id, 
                      predictPositions=seq(0,1500,by=10))

## ----eval=FALSE----------------------------------------------------------
#  help(Bchronology)

## ----eval=FALSE----------------------------------------------------------
#  summary(GlenOut)

## ------------------------------------------------------------------------
summary(GlenOut, type='convergence')
summary(GlenOut, type='outliers')

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
plot(GlenOut,
     main="Glendalough",
     xlab='Age (cal years BP)',
     ylab='Depth (cm)',
     las=1)

## ----results='hide'------------------------------------------------------
predictAges = predict(GlenOut, 
                      newPositions = c(150,725,1500), 
                      newPositionThicknesses=c(5,0,20))

## ---- results ='hide'----------------------------------------------------
acc_rate = summary(GlenOut, type = 'acc_rate')

## ---- eval=FALSE---------------------------------------------------------
#  plot(acc_rate[,'age_grid'], acc_rate[,'50%'], type='l', ylab = 'cm per year', xlab = 'Age (k cal years BP)', ylim = range(acc_rate[,-1]))
#  lines(acc_rate[,'age_grid'], acc_rate[,'2.5%'], lty='dotted')
#  lines(acc_rate[,'age_grid'], acc_rate[,'97.5%'], lty='dotted')

## ---- eval=FALSE---------------------------------------------------------
#  sed_rate = summary(GlenOut, type = 'sed_rate', useExisting = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  plot(sed_rate[,'position_grid'], sed_rate[,'50%'], type='l', ylab = 'Years per cm', xlab = 'Depth (cm)', ylim = range(sed_rate[,-1]))
#  lines(sed_rate[,'position_grid'], sed_rate[,'2.5%'], lty='dotted')
#  lines(sed_rate[,'position_grid'], sed_rate[,'97.5%'], lty='dotted')

## ------------------------------------------------------------------------
data(TestChronData)
data(TestRSLData)

## ----messages=FALSE, results='hide', eval=FALSE--------------------------
#  RSLrun = Bchronology(ages=TestChronData$ages,
#                       ageSds=TestChronData$ageSds,
#                       positions=TestChronData$position,
#                       positionThicknesses=TestChronData$thickness,
#                       ids=TestChronData$id,
#                       calCurves=TestChronData$calCurves,
#                       predictPositions=TestRSLData$Depth)
#  RSLrun2 = BchronRSL(RSLrun,
#                      RSLmean=TestRSLData$RSL,
#                      RSLsd=TestRSLData$Sigma,
#                      degree=3)

## ---- eval=FALSE---------------------------------------------------------
#  summary(RSLrun2, type = 'RSL', age_grid = seq(0, 2000, by  = 250))
#  plot(RSLrun2, type = 'RSL', main = 'Relative sea level plot')
#  plot(RSLrun2, type = 'rate', main = 'Rate of RSL change')

## ----results='hide', eval=FALSE------------------------------------------
#  data(Sluggan)
#  SlugDens = BchronDensity(ages=Sluggan$ages,
#                           ageSds=Sluggan$ageSds,
#                           calCurves=Sluggan$calCurves,
#                           numMix=50)

## ---- eval=FALSE---------------------------------------------------------
#  plot(SlugDens,xlab='Age (cal years BP)')

## ----eval=FALSE----------------------------------------------------------
#  SlugDensFast = BchronDensityFast(ages=Sluggan$ages,
#                                   ageSds=Sluggan$ageSds,
#                                   calCurves=Sluggan$calCurves)

## ------------------------------------------------------------------------
# Load in the calibration curve with:
intcal09 = read.table('http://www.radiocarbon.org/IntCal09%20files/intcal09.14c',sep=',')
# Run CreateCalCurve
CreateCalCurve(name='intcal09',cal_ages=intcal09[,1],uncal_ages=intcal09[,2],one_sigma=intcal09[,3])

## ------------------------------------------------------------------------
age_09 = BchronCalibrate(age=15500,ageSds=150,calCurves = 'intcal09',ids='My Date')
age_13 = BchronCalibrate(age=15500,ageSds=150,calCurves = 'intcal13')
plot(age_09)
lines(age_13$Date1$ageGrid,age_13$Date1$densities,col='red')
legend('topleft',legend=c('intcal09','intcal13'),col=c('black','red'),lty=1)

