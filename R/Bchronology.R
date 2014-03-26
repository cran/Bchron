# Function to create Bchron chronologies
Bchronology = function(ages,ageSds,positions,positionThicknesses=rep(0,length(ages)),calCurves=rep('intcal13',length(ages)),ids=NULL,outlierProbs=rep(0.01,length(ages)),predictPositions=seq(min(positions),max(positions),length=100),pathToCalCurves=system.file('data',package='Bchron'),iterations=10000,burn=2000,thin=8,extractDate=1950-as.numeric(format(Sys.time(),"%Y")),maxExtrap=500,thetaMhSd=0.5,muMhSd=0.1,psiMhSd=0.1,ageScaleVal=1000,positionScaleVal=100) {
  
# Notation:
# theta are the calibrated ages of ages 1 to n (not necessarily radiocarbon)
# phi are the outlier indicators (1=TRUE or 0=FALSE) for date i
# mu,psi are the Compound Poisson-Gamma parameters controlling sedimentation

# Check positions don't overlap
n = length(ages)
depthLow = positions-0.5*positionThicknesses
depthHigh = positions+0.5*positionThicknesses
for(i in 2:n) if(depthLow[i]-depthHigh[i-1]<0 & all(positionThicknesses==0)) stop('Depth layers identical with no thickness errors - not supported in Bchron. Check thickness errors')

# Re-scale the predicted positions
predictPositionsRescaled = predictPositions/positionScaleVal

# Check positions are in order
o = order(positions)
if(any(positions[o]!=positions)) {
  warning("positions not given in order - re-ordering")
  ages = ages[o]
  ageSds = ageSds[o]
  positions=positions[o]
  positionThicknesses=positionThicknesses[o]
  ids = ids[o]
  outlierProbs = outlierProbs[o]
}

# First thing to do is to calibrate the ages with the two different degrees of freedom
dfs = c(1,100) # corresponding to phi equals 1 and 0 respectively
x.df1 = BchronCalibrate(ages=ages,ageSds=ageSds,calCurves=calCurves,positions=positions,ids=ids,pathToCalCurves=pathToCalCurves,eps=0,dfs=rep(dfs[1],length(ages)))
x.df2 = BchronCalibrate(ages=ages,ageSds=ageSds,calCurves=calCurves,positions=positions,ids=ids,pathToCalCurves=pathToCalCurves,eps=0,dfs=rep(dfs[2],length(ages)))

# Get current positions and their order
currPositions = sort(jitter(positions/positionScaleVal))
diffPosition = diff(currPositions)
do = order(currPositions)

# For any calibration curves that don't start at 0, we need an offset to enable fast lookup - only supported one so far is normal
offset=rep(0,length=n)
for(i in 1:n) {
   offset[i] = ifelse(x.df1[[i]]$calCurve == 'normal',61,0)
}

# Starting values
theta = vector(length=n)
# Make sure no theta values are identical 
badThetas = TRUE
while(badThetas) {
  for(j in 1:n) theta[j] = round(rnorm(1,x.df2[[j]]$ageGrid[match(max(x.df2[[j]]$densities),x.df2[[j]]$densities)]/ageScaleVal,sd=ageSds[j]/ageScaleVal),3)
  theta=sort(theta)
  if(all(diff(theta)!=0)) badThetas= FALSE
}
phi = rep(0,length(theta))
p = 1.2
mu = abs(rnorm(1,mean=mean(diff(theta))/mean(diffPosition),sd=1))
psi = abs(rnorm(1,1,1))

# Tranformed values (used for interpolation)
alpha = (2-p)/(p-1)

# Storage
remaining=(iterations-burn)/thin
thetaStore = phiStore = matrix(ncol=length(theta),nrow=remaining)
muStore = psiStore = vector(length=remaining)
thetaPredict = matrix(ncol=length(predictPositions),nrow=remaining)

# Some C functions which are useful 
#################################################

rtruncn.sig = signature(a='numeric',b='numeric',x='numeric')
rtruncncode = '
double A, B;
double maxA, maxB, maxR, r2, r, th, u, v, accept=0.0;
A = atan(*a);
B = atan(*b);
maxA = exp(-pow(*a,2)/4)/cos(A);
maxB = exp(-pow(*b,2)/4)/cos(B);
maxR = fmax2(maxA, maxB);
if((*a<1) && (*b>-1)) maxR = exp(-0.25)*sqrt(2.0);
while (accept==0) {
r2 = runif(0.0,1.0);
r = sqrt(r2)*maxR;
th = runif(A,B);
u = r*cos(th);
*x = tan(th);
accept = ((pow(*x,2)) < (log(u)*-4));
}
'

truncatedWalk.sig = signature(old='numeric',sd='numeric',low='numeric',high='numeric',newvalue='numeric')
truncatedWalkcode = '
double lowlimold, upplimold, y;
lowlimold = (*low - *old)/ *sd;
upplimold = (*high - *old)/ *sd;
rtruncn(&lowlimold, &upplimold,&y);
*newvalue = *old + *sd*y;
'

truncatedrat.sig = signature(old='numeric',sd='numeric',low='numeric',high='numeric',newvalue='numeric',ratio='numeric')
truncatedratcode = '
double lowlimold, upplimold, lowlimnew, upplimnew, plowold, puppold, plownew, puppnew;
lowlimold = (*low - *old)/ *sd;
upplimold = (*high - *old)/ *sd;
lowlimnew = (*low - *newvalue)/ *sd;
upplimnew = (*high - *newvalue)/ *sd;
plowold = pnorm(lowlimold,0.0,1.0,1,0);
puppold = pnorm(upplimold,0.0,1.0,1,0);
plownew = pnorm(lowlimnew,0.0,1.0,1,0);
puppnew = pnorm(upplimnew,0.0,1.0,1,0);
*ratio = (puppold - plowold)/(puppnew - plownew);
'

truncatedWalk.Funs = cfunction(list(rtruncn=rtruncn.sig,truncatedWalk=truncatedWalk.sig,truncatedrat=truncatedrat.sig),list(rtruncncode,truncatedWalkcode,truncatedratcode),language='C',convention='.C',verbose=FALSE,includes='#include <Rmath.h>')
truncatedWalk.raw = truncatedWalk.Funs[['truncatedWalk']]
truncatedrat.raw = truncatedWalk.Funs[['truncatedrat']]

truncatedWalk = function(old,sd,low,high) {
  new=truncatedWalk.raw(old,sd,low,high,0)$newvalue
  rat=truncatedrat.raw(old,sd,low,high,new,0)$ratio
  if(is.nan(rat)) rat=1 # Useful for when proposed value and new value are identical
  return(list(new=new,rat=rat))
}

# C functions for tweedie
dtweedielogwsmallp.sig = signature(y='numeric',phi='numeric',power='numeric',logw='numeric')
dtweedielogwsmallp.code = '
  double p,a,a1,r,drop=37,logz,jmax,j,cc,wmax,estlogw,oldestlogw;
  int hij,lowj;
  
  if (*power < 1) error("Error - power<1!");
  if (*power > 2) error("Error - power>2!");
  if (*phi <= 0) error("Error - phi<=0!");
  if (*y <= 0) error("Error - y<=0!");
  p = *power;
  a = (2 - p)/(1 - p);
  a1 = 1 - a;
  r = -a * log(*y) + a * log(p - 1) - a1 * log(*phi) - log(2 - p);
  logz = r;
  
  jmax = (pow(*y,(2 - p)))/(*phi * (2 - p));
  j = fmax2(1, jmax);
  cc = logz + a1 + a * log(-a);
  wmax = a1 * jmax;
  estlogw = wmax;
  while (estlogw > (wmax - drop)) 
{
  j = j + 2;
  estlogw = j * (cc - a1 * log(j));
}
  
  hij = (int)ceil(j);
  logz = r;
  jmax = pow(*y,(2 - *power))/(*phi * (2 - *power));
  j = fmax2(1, jmax);
  wmax = a1 * jmax;
  estlogw = wmax;
  while ((estlogw > (wmax - drop)) && (j >= 2)) 
{
  j = fmax2(1, j - 2);
  oldestlogw = estlogw;
  estlogw = j * (cc - a1 * log(j));
}
  lowj = (int)fmax2(1, floor(j));
  
  double newj[hij-lowj+1];
  int k;
  for(k=0;k<(hij-lowj+1);k++) newj[k] = lowj+k;
  
  double g[hij-lowj+1]; 
  for(k=0;k<hij-lowj+1;k++) g[k] = lgamma(newj[k]+1)+lgamma(-a*newj[k]);
  
  double A[hij-lowj+1];
  for(k=0;k<hij-lowj+1;k++) A[k] = r*(double)newj[k]-g[k];
  
  double m=fmax2(A[0],hij-lowj+1);
  for(k=0;k<(hij-lowj+1);k++) m = fmax2(A[k],hij-lowj+1);

  double we[hij-lowj+1];
  for(k=0;k<hij-lowj+1;k++) we[k] = exp(A[k]-m);
  double sumwe=0;
  for(k=0;k<hij-lowj+1;k++) sumwe+=we[k];
  *logw=log(sumwe)+m;
'

dtweedieseriessmallp.sig = signature(power='numeric',y='numeric',mu='numeric',phi='numeric',f='numeric')
dtweedieseriessmallp.code = 
'  
  double logw;
  dtweedielogwsmallp(y,phi,power,&logw);
  double tau = *phi*(*power-1)*pow(*mu,*power-1);
  double lambda = pow(*mu,2-*power)/(*phi*(2-*power));
  double logf = -*y/tau-lambda-log(*y)+logw;
  *f = exp(logf);
'

dtweediep1.sig = signature(y='numeric',power='numeric',mu='numeric',phi='numeric',fTplus='numeric')
dtweediep1.code = 
'
  // Calculates the density of a tweedie plus one random variable
  double eps = 0.00000001;
  double lambda2 = pow(*mu,2-*power)/(*phi*(2-*power))-eps;
  double alpha = (2-*power)/(*power-1);
  double beta = 1/(*phi*(*power-1)*pow(*mu,*power-1));
  
  double mu2 = alpha*lambda2/beta;
  double phi2 = (alpha+1)/(pow(lambda2*alpha,(1/(alpha+1)))*pow(beta,(alpha/(alpha+1))));
  
  double fTplus1,fTplus2,fTplus3;
  dtweedieseriessmallp(power,y,mu,phi,&fTplus1);
  dtweedieseriessmallp(power,y,mu,phi,&fTplus2);
  dtweedieseriessmallp(power,y,&mu2,&phi2,&fTplus3);
  *fTplus = fTplus1+(1/eps)*(fTplus2-fTplus3);
'

Tweedie.Funs = cfunction(list(dtweedielogwsmallp=dtweedielogwsmallp.sig,dtweedieseriessmallp=dtweedieseriessmallp.sig,dtweediep1=dtweediep1.sig),list(dtweedielogwsmallp.code,dtweedieseriessmallp.code,dtweediep1.code),language='C',convention='.C',verbose=FALSE,includes='#include <Rmath.h>')

dtweediep1.raw = Tweedie.Funs[['dtweediep1']]
dtweediep1 = Vectorize(function(y,p,mu,phi) {
  return(dtweediep1.raw(y,p,mu,phi,0)$fTplus)
})

# Some C functions to do prediction and interpolation
predictInterp.sig = signature(alpha='numeric',lambda='numeric',beta='numeric',predictPositions='numeric',NpredictPositions='integer',diffPositionj='numeric',currPositionsj='numeric',currPositionsjp1='numeric',thetaj='numeric',thetajp1='numeric',predvals='numeric')
predictInterp.code = '
// Runs the prediction code when we are interpolating between two positions
int Nd = rpois((*lambda)*(*diffPositionj));
int i;
double depthEvents[Nd];
for(i=0;i<Nd;i++) depthEvents[i] = runif(*currPositionsj,*currPositionsjp1);
R_rsort(depthEvents,Nd);
double timeEventsUnsc[Nd+1],timeEventsSum=0.0;
for(i=0;i<Nd+1;i++) timeEventsUnsc[i] = rgamma(*alpha,1/(*beta));
for(i=0;i<Nd+1;i++) timeEventsSum += timeEventsUnsc[i];
double timeEvents[Nd+1];
for(i=0;i<Nd+1;i++) timeEvents[i] = (*thetajp1-*thetaj)*timeEventsUnsc[i]/timeEventsSum;
double timeEventsCumsum[Nd+1],allTimeEvents[Nd+2];
timeEventsCumsum[0] = 0.0;
for(i=1;i<Nd+1;i++) timeEventsCumsum[i] = timeEventsCumsum[i-1] + timeEvents[i];
for(i=0;i<Nd+1;i++) allTimeEvents[i] = timeEventsCumsum[i]+*thetaj;
allTimeEvents[Nd+1] = *thetajp1;
double allDepthEvents[Nd+2];
allDepthEvents[0] = *currPositionsj;
for(i=1;i<Nd+1;i++) allDepthEvents[i] = depthEvents[i-1];
allDepthEvents[Nd+1] = *currPositionsjp1;

int Ndp2 = Nd+2;
for(i=0;i<*NpredictPositions;i++) {
  linInterp(&Ndp2,&predictPositions[i],allDepthEvents,allTimeEvents,&predvals[i]);
}
'

linInterp.sig = signature(n='integer',newx='numeric',x='numeric',y='numeric',ans='numeric')
linInterp.code = '
// Try to predict the value newx from the vectors x and y
for(int i=0; i<*n-1; i++) {
  if(((*newx >= x[i]) & (*newx <= x[i+1])) | ((*newx <= x[i]) & (*newx >= x[i+1]))) {
    *ans = y[i] + ((*newx-x[i])/(x[i+1]-x[i]))*(y[i+1]-y[i]);
    if(*newx==x[i]) *ans = y[i];
  }        
}
'

predictExtrapUp.sig = signature(alpha='numeric',lambda='numeric',beta='numeric',predictPositions='numeric',NpredictPositions='integer',currPositions1='numeric',theta1='numeric',maxExtrap='integer',extractDate='numeric',predvals='numeric')
predictExtrapUp.code = '
// Runs the prediction code when we are extrapolating up beyond the first date
int bad=1,count=0,i;
double depthEvents[*maxExtrap],timeEvents[*maxExtrap];
depthEvents[0] = *currPositions1;
timeEvents[0] = *theta1;
while(bad==1) {
  for(i=1;i<*maxExtrap;i++) {
    depthEvents[i] =  depthEvents[i-1]-rexp(1/(*lambda));
    timeEvents[i] =  timeEvents[i-1]-rgamma(*alpha,1/(*beta));
  }
  for(i=0;i<*NpredictPositions;i++) {
    linInterp(maxExtrap,&predictPositions[i],depthEvents,timeEvents,&predvals[i]);
  }
  count+=1;
  bad=0;
  for(i=0;i<*NpredictPositions;i++) {
    if(predvals[i]<*extractDate) bad=1;
  }
  if(count==50) {
    for(i=0;i<*NpredictPositions;i++) {
      if(predvals[i]<*extractDate) predvals[i] = *extractDate;
    }    
    bad=0;
    warning("Unable to find suitable chronologies for top of core - truncated to date of extraction");
  }
}
'     

predictExtrapDown.sig = signature(alpha='numeric',lambda='numeric',beta='numeric',predictPositions='numeric',NpredictPositions='integer',currPositionsn='numeric',thetan='numeric',maxExtrap='integer',predvals='numeric')
predictExtrapDown.code = '
// Runs the prediction code when we are extrapolating down below the bottom date
double depthEvents[*maxExtrap],timeEvents[*maxExtrap];
int i;
depthEvents[0] = *currPositionsn;
timeEvents[0] = *thetan;
for(i=1;i<*maxExtrap;i++) {
  depthEvents[i] =  depthEvents[i-1]+rexp(1/(*lambda));
  timeEvents[i] =  timeEvents[i-1]+rgamma(*alpha,1/(*beta));
}
for(i=0;i<*NpredictPositions;i++) {
  linInterp(maxExtrap,&predictPositions[i],depthEvents,timeEvents,&predvals[i]);
}
'   


predict.Funs = cfunction(list(linInterp=linInterp.sig,predictInterp=predictInterp.sig,predictExtrapUp=predictExtrapUp.sig,predictExtrapDown=predictExtrapDown.sig),list(linInterp.code,predictInterp.code,predictExtrapUp.code,predictExtrapDown.code),language='C',convention='.C',verbose=FALSE,includes=c('#include <Rmath.h>','#include <R.h>'))

predictInterp.raw = predict.Funs[['predictInterp']]
predictInterp = function(alpha,lambda,beta,predictPositions,diffPositionj,currPositionsj,currPositionsjp1,thetaj,thetajp1) {
  return(predictInterp.raw(alpha,lambda,beta,predictPositions,length(predictPositions),diffPositionj,currPositionsj,currPositionsjp1,thetaj,thetajp1,rep(0,length(predictPositions)))$predvals)
}
predictExtrapUp.raw = predict.Funs[['predictExtrapUp']]
predictExtrapUp = function(alpha,lambda,beta,predictPositions,currPositions1,theta1,maxExtrap,extractDate) {
  return(predictExtrapUp.raw(alpha,lambda,beta,predictPositions,length(predictPositions),currPositions1,theta1,maxExtrap,extractDate,rep(0,length(predictPositions)))$predvals)
}
predictExtrapDown.raw = predict.Funs[['predictExtrapDown']]
predictExtrapDown = function(alpha,lambda,beta,predictPositions,currPositionsn,thetan,maxExtrap) {
  return(predictExtrapDown.raw(alpha,lambda,beta,predictPositions,length(predictPositions),currPositionsn,thetan,maxExtrap,rep(0,length(predictPositions)))$predvals)
}

# End of C functions

# Main iteration loop
#################################################

pb = utils::txtProgressBar(min = 1, max = iterations, style = 3,width=60,title='Running Bchronology...')
for(i in 1:iterations) {
  utils::setTxtProgressBar(pb, i)
  
  if(any(positionThicknesses>0) & i>0.5*burn & i%%thin==0) {
    currPositions = runif(n,positions/positionScaleVal-0.5*positionThicknesses/positionScaleVal,positions/positionScaleVal+0.5*positionThicknesses/positionScaleVal)
    # Get date order so I can preserve things if they change around
    do = order(currPositions)
    diffPosition = diff(currPositions[do])
    theta[do]=sort(theta)
  }
  
  # If we're in the right place put things in storage
  if(i>burn & i%%thin==0) {
    ind = (i-burn)/thin
    thetaStore[ind,] = theta*ageScaleVal
    phiStore[ind,] = phi
    muStore[ind] = mu
    psiStore[ind] = psi
        
    # Run interpolation/extrapolation stage
    lambda = (mu^(2-p))/(psi*(2-p))
    beta = 1/(psi*(p-1)*(mu^(p-1)))
    
    # First interpolation
    for(j in 1:(n-1)) {
      # Find which positions we need to interpolate for
      depthIndRange = which(predictPositionsRescaled>=currPositions[do[j]] & predictPositionsRescaled<=currPositions[do[j+1]])
      thetaPredict[ind,depthIndRange] = round(predictInterp(alpha,lambda,beta,predictPositionsRescaled[depthIndRange],diffPosition[j],currPositions[do[j]],currPositions[do[j+1]],theta[do[j]],theta[do[j+1]]),3)
    }
  
    # Extrapolate up to to top depth
    if(any(predictPositionsRescaled<currPositions[1])) {
      depthIndRange = which(predictPositionsRescaled<=currPositions[1])
      thetaPredict[ind,depthIndRange] = round(predictExtrapUp(alpha,lambda,beta,predictPositionsRescaled[depthIndRange],currPositions[1],theta[1],maxExtrap,extractDate/ageScaleVal),3)
    }
    
    # Extrapolate below bottom depth
    if(any(predictPositionsRescaled>=currPositions[n])) {
      depthIndRange = which(predictPositionsRescaled>=currPositions[n])
      thetaPredict[ind,depthIndRange] = round(predictExtrapDown(alpha,lambda,beta,predictPositionsRescaled[depthIndRange],currPositions[n],theta[n],maxExtrap),3)
    }

    if(any(is.na(thetaPredict[ind,]))) stop("Errors in predicted ages. Check you are not extrapolating too far away from dated levels. If you must run this core with these predicted ages, set maxExtrap to a larger value (e.g. 1000)")
  }

  # Update theta
  for(j in 1:n) {
    thetaNewAll = truncatedWalk(theta[do[j]],thetaMhSd,ifelse(j==1,ifelse(x.df1[[j]]$calCurve=='normal',extractDate/ageScaleVal,0),theta[do[j-1]]+0.001),ifelse(j==n,100000,theta[do[j+1]]-0.001))
    thetaNew = round(thetaNewAll$new,3)
    # Calculate ratio
    if(phi[do[j]]==0) {
      currDens = x.df2[[do[j]]]$densities
    } else {
      currDens = x.df1[[do[j]]]$densities
    }
    thetaNewMatch = as.integer(thetaNew*ageScaleVal+offset[do[j]])+1      
    thetaNewLogDens = max(log(currDens[thetaNewMatch]),-1000000)
    priorNewLogDens = ifelse(j==1,0,log(dtweediep1(thetaNew-theta[do[j-1]],p,mu*diffPosition[j-1],psi/(diffPosition[j-1]^(p-1)))))+ifelse(j==n,0,log(dtweediep1(theta[do[j+1]]-thetaNew,p,mu*(diffPosition[j]),psi/(diffPosition[j])^(p-1))))
    thetaMatch = as.integer(theta[do[j]]*ageScaleVal+offset[do[j]])+1
    thetaLogDens = max(log(currDens[thetaMatch]),-1000000)
    priorLogDens = ifelse(j==1,0,log(dtweediep1(theta[do[j]]-theta[do[j-1]],p,mu*(diffPosition[j-1]),psi/(diffPosition[j-1])^(p-1))))+ifelse(j==n,0,log(dtweediep1(theta[do[j+1]]-theta[do[j]],p,mu*(diffPosition[j]),psi/(diffPosition[j])^(p-1))))
    
    logRtheta = thetaNewLogDens - thetaLogDens + priorNewLogDens - priorLogDens + log(thetaNewAll$rat)
    if(runif(1)<exp(logRtheta)) theta[do[j]] = thetaNew      
  }
  
  # Update phi
  for(j in 1:n) {
    phiNew = sample(0:1,1)
    if(phiNew != phi[do[j]]) {
      if(phi[do[j]]==0) {
        currDens = x.df2[[do[j]]]$densities
      } else {
        currDens = x.df1[[do[j]]]$densities
      }
      thetaMatch = as.integer(theta[do[j]]*ageScaleVal+offset[j])+1
      thetaLogDens = max(log(currDens[thetaMatch]),-1000000)
      if(phiNew==0) {
        newDens = x.df2[[do[j]]]$densities
      } else {
        newDens = x.df1[[do[j]]]$densities
      }
      thetaMatch = as.integer(theta[do[j]]*ageScaleVal+offset[j])+1
      thetaNewLogDens = max(log(newDens[thetaMatch]),-1000000)
      
      logRphi = thetaNewLogDens- thetaLogDens + dbinom(phiNew,1,outlierProbs[do[j]],log=TRUE) - dbinom(phi[do[j]],1,outlierProbs[do[j]],log=TRUE)
      
      if(runif(1)<exp(logRphi)) phi[do[j]] = phiNew
    }
  }
  
  # Update mu
  muNewAll = truncatedWalk(mu,muMhSd,0,1e5)
  muNew = muNewAll$new

  logRmu = sum(log(dtweediep1(diff(theta[do]),p,muNew*diffPosition,psi/(diffPosition)^(p-1)))) - sum(log(dtweediep1(diff(theta[do]),p,mu*diffPosition,psi/(diffPosition)^(p-1)))) + log(muNewAll$rat)
  if(runif(1)<exp(logRmu)) mu = muNew
  
  # Update psi
  psiNewAll = truncatedWalk(psi,psiMhSd,0,1e5)
  psiNew = psiNewAll$new
  
  logRpsi = sum(log(dtweediep1(diff(theta[do]),p,mu*diffPosition,psiNew/(diffPosition)^(p-1)))) - sum(log(dtweediep1(diff(theta[do]),p,mu*diffPosition,psi/(diffPosition)^(p-1)))) + log(psiNewAll$rat)
  if(runif(1)<exp(logRpsi)) psi = psiNew
  
}

# Return everything
out = list(theta=thetaStore,phi=phiStore,mu=muStore,psi=psiStore,thetaPredict=ageScaleVal*thetaPredict,predictPositions=predictPositions,calAges=x.df2)
class(out) = 'BchronologyRun'
return(out)
  
}