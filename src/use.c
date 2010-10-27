// A couple of files that more than one thingy uses
#include<R.h>
#include<Rmath.h>
#include<stdlib.h>

int diff(double *arr,int *len,double *retarr)
{
// this function takes a one-dimensional array arr and its length len, and returns the differenced
// vector retarr of length len-1

int i;
for(i=0;i<*len-1;i++)
{
	retarr[i] = arr[i+1]-arr[i];
}
return(0);

}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

////////////////////// Some other functions //////////////////////////////////


int GetLengthCurrentDepths(double depthlow,double depthhigh,double *ddepth,int length)
{
// This function finds the sizes of the array required to store the current depths

// Need to loop through the depths to get the number of things in ddepths which are
// between depthlow and depthhigh

int i,count=0;
for(i=0;i<length;i++)
{
    if((ddepth[i]<=depthhigh) & (ddepth[i]>=depthlow)) count = count+1;

}

return(count);
}

int GetCurrentDepths(double depthlow,double depthhigh,double *ddepth,int length,double *current)
{
// This function takes currentdepths and fills it with all the ddepths between depthlow and
// depthhigh

int i,count=0;
for(i=0;i<length;i++)
{
    if((ddepth[i]<=depthhigh) & (ddepth[i]>=depthlow))
    {
      current[count] = ddepth[i];
      count = count+1;  
    }
    
}

return(0);
}

int GetCurrentDepthRows(double depthlow,double depthhigh,double *ddepth,int length,int *rows)
{
// This function takes currentdepths and fills it with all the ddepths between depthlow and
// depthhigh

int i,count=0;
for(i=0;i<length;i++)
{
    if((ddepth[i]<=depthhigh) & (ddepth[i]>=depthlow))
    {
      rows[count] = i;
      count = count+1;  
    }
    
}

return(0);
}

double linearinterp(int n, double newx, double *a, double *b)
{
    double newvalue;
    int i;

//condition is where y lies between the two closests approximations to 
//it in the cal curve
   
    for(i=0; i<n-1; i++)
    {
        if (((newx >= a[i]) & (newx <= a[i+1])) | ((newx <= a[i]) & (newx >= a[i+1])))
        {
                newvalue = b[i] + ((newx-a[i])/(a[i+1]-a[i]))*(b[i+1]-b[i]);
                if(newx==a[i]) newvalue = b[i];
                return(newvalue);
                //break;
        }        
    }
  
  return(-999.0);
}

double densinterp(int nrows, int ncols, double newx, double a[nrows][ncols],double b[nrows][ncols],int choosecol,double outval)
{
    double newvalue;
    int i;

	// This function is similar to the linear interp function but instead forces values
	// outside the range to be the value outval (usually zero).

	if( (newx < a[0][choosecol]) | (newx > a[nrows-1][choosecol]) ) {
 		return(outval);
	} else {
	    for(i=0; i<nrows-1; i++)
    	{
        	if (((newx >= a[i][choosecol]) & (newx <= a[i+1][choosecol])) | ((newx <= a[i][choosecol]) & (newx >= a[i+1][choosecol])))
	        {
    	            newvalue = b[i][choosecol] + ((newx-a[i][choosecol])/(a[i+1][choosecol]-a[i][choosecol]))*(b[i+1][choosecol]-b[i][choosecol]);
        	        //if(newx==a[i]) newvalue = b[i];
            	    return(newvalue);
                	//break;
	        }        
    	}
	}
  
  Rprintf("Weird density interpolation. Don't trust these results!\n");
  return(-999.0);
}

double Max2 (double a, double b)
{
	// find the max of 2 numbers
   double larger;
   if (a > b)
      larger = a;
   else
      larger = b;
   return larger;
}

//rtruncn function:
double rtruncn (double a, double b)
{
    double A, B;
    double maxA, maxB, maxR, r2, r, th, u, v, x, accept=1.0;
    
    A = atan(a);
    B = atan(b);
    
    maxA = exp(-pow(a,2)/4)/cos(A);
    maxB = exp(-pow(b,2)/4)/cos(B);
    maxR = Max2(maxA, maxB);

    if((a<1) && (b>-1)) maxR = exp(-0.25)*sqrt(2.0);

    while (accept!=0)
    {
        r2 = runif(0.0,1.0);
        r = sqrt(r2)*maxR;
        th = runif(A,B);
        u = r*cos(th);
        v = r*sin(th);
        x = tan(th);
        accept = ((pow(x,2)) < (log(u)*-4));
    }
    return x;

}        

//truncated normal function:
double truncatedwalk (double old, double sd, double low, double high)
{
    double lowlimold, upplimold, y, newvalue;
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    y = rtruncn(lowlimold, upplimold);
    newvalue = old + sd*y;
           
    return newvalue;
}

//truncated normal ratio function:
double truncatedrat (double old, double sd, double low, double high, double newvalue)
{
    double lowlimold, upplimold, lowlimnew, upplimnew, plowold, puppold, plownew, puppnew, ratio;
    
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    lowlimnew = (low - newvalue)/sd;
    upplimnew = (high - newvalue)/sd;
    plowold = pnorm(lowlimold,0.0,1.0,1,0);
    puppold = pnorm(upplimold,0.0,1.0,1,0);
    plownew = pnorm(lowlimnew,0.0,1.0,1,0);
    puppnew = pnorm(upplimnew,0.0,1.0,1,0);
    ratio = (puppold - plowold)/(puppnew - plownew);
    return ratio;        
}

double Max(double *Numbers, int Count)
{
	// Find the maximum of a sequence of numbers
	double Maximum;
	Maximum = Numbers[0];

	for(int i = 0; i < Count; i++)
		if( Maximum < Numbers[i] )
			Maximum = Numbers[i];

	return Maximum;
}

double Min(double *Numbers, int Count)
{
	// Find the maximum of a sequence of numbers
	double Minimum;
	Minimum = Numbers[0];

	for(int i = 0; i < Count; i++)
		if( Minimum > Numbers[i] )
			Minimum = Numbers[i];

	return Minimum;
}

int seq(double from,double to,double len,double *sequence)
{
	// Create a sequence of numbers from 'from' to 'to' of length 'len'
	// Simple huh?
	
	double by = (to-from)/(len-1);
	int i;
	for(i=0;i<len;i++) 
		sequence[i] = from + i*by;

	return(0);

}

double dtweedielogwsmallp(double y, double phi, double power)
{
	// Matches the R function of the same name. Oh, actually it's called 

	double p,a,a1,r,drop=37,logz,jmax,j,cc,wmax,estlogw,oldestlogw;
	int hij,lowj;

    if (power < 1) 
        exit(-99);
	if (power > 2) 
		exit(-99);
    if (phi <= 0)
		exit(-99);
    if (y <= 0)
		exit(-99);
    p = power;
    a = (2 - p)/(1 - p);
    a1 = 1 - a;
    r = -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p);
    logz = r;
	
    jmax = (pow(y,(2 - p)))/(phi * (2 - p));
    j = Max2(1, jmax);
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
    jmax = pow(y,(2 - power))/(phi * (2 - power));
    j = Max2(1, jmax);
    wmax = a1 * jmax;
    estlogw = wmax;
    while ((estlogw > (wmax - drop)) && (j >= 2)) 
	{
        j = Max2(1, j - 2);
        oldestlogw = estlogw;
        estlogw = j * (cc - a1 * log(j));
    }
    lowj = (int)Max2(1, floor(j));

	double newj[hij-lowj+1];
    seq(lowj, hij,(hij-lowj+1),newj);
    
	double g[hij-lowj+1]; 
	int k;
	for(k=0;k<hij-lowj+1;k++) g[k] = lgamma(newj[k]+1)+lgamma(-a*newj[k]);
	
	double A[hij-lowj+1];
	for(k=0;k<hij-lowj+1;k++) A[k] = r*(double)newj[k]-g[k];
	
	double m=Max(A,hij-lowj+1);
    double we[hij-lowj+1];
	for(k=0;k<hij-lowj+1;k++) we[k] = exp(A[k]-m);
	double sumwe=0;
	for(k=0;k<hij-lowj+1;k++) sumwe+=we[k];
	double logw=log(sumwe)+m;

	return(logw);

}

double dtweedieseriessmallp(double power,double y, double mu, double phi)
{

// This function matches the R function of the same name (with a few dots in it though)

double logw = dtweedielogwsmallp(y,phi,power);
double tau = phi*(power-1)*pow(mu,power-1);
double lambda = pow(mu,2-power)/(phi*(2-power));
double logf = -y/tau-lambda-log(y)+logw;
double f = exp(logf);

return(f);

}

double dtweediep1(double y, double power, double mu, double phi)
{
// Same as my R function
// Calculates the density of a tweedie plus one random variable

double eps = 0.00000001;
double lambda2 = pow(mu,2-power)/(phi*(2-power))-eps;
double alpha = (2-power)/(power-1);
double beta = 1/(phi*(power-1)*pow(mu,power-1));

double mu2 = alpha*lambda2/beta;
double phi2 = (alpha+1)/(pow(lambda2*alpha,(1/(alpha+1)))*pow(beta,(alpha/(alpha+1))));

double fTplus = dtweedieseriessmallp(power,y,mu,phi)
	+(1/eps)*(dtweedieseriessmallp(power,y,mu,phi)
	-dtweedieseriessmallp(power,y,mu2,phi2));

return(fTplus);

}

double UpdateMCMC(double newloglik,double oldloglik,double newval,double oldval,double rat)
{
// Function to update MCMC when given log likelihoods

double u,mh;
u = runif(0.0,1.0);
mh = exp(newloglik - oldloglik)*rat;
if (u < mh)
{
	return(newval);
} else 
{
	return(oldval);
}

}

int fact(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * fact(number - 1);
	return temp;
}

double dlogbinom(int x,int n, double p)
{
double dlogbinom = log(fact(n))-log(fact(x))-log(fact(n-x))+x*log(p)+(n-x)*log(1-p);
return dlogbinom;

}

int compare_doubles (const double *a, const double *b)
{
  double temp = *a - *b;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

void ReOrder(double *currentdepths,double *cage,double *sd,double *depth,double *thick, double *outprob1, double *outprob2, int *type, int ndets)
{
// This function takes the current depths when they are not in order and places
// everything back in the right order

int Order[ndets];
double sorteddepths[ndets];
double currentdepthst[ndets],caget[ndets],sdt[ndets],thickt[ndets],outprob1t[ndets];
double deptht[ndets],outprob2t[ndets],typet[ndets];
int i;

// First create sorted depths
for(i=0;i<ndets;i++) {
    sorteddepths[i] = currentdepths[i];
    Order[i]=i;
}

// Now sort
//R_rsort(sorteddepths,ndets);
R_qsort_I(sorteddepths,Order,1,ndets);

for(i=0;i<ndets;i++) { 
    currentdepthst[i] = currentdepths[Order[i]];
    caget[i] = cage[Order[i]];
    sdt[i] = sd[Order[i]];
    deptht[i] = depth[Order[i]];
    thickt[i] = thick[Order[i]];
    outprob1t[i] = outprob1[Order[i]];
    outprob2t[i] = outprob2[Order[i]];
    typet[i] = type[Order[i]];
}

for(i=0;i<ndets;i++) { 
    currentdepths[i] = currentdepthst[i];
    cage[i] = caget[i];
    sd[i] = sdt[i];
    depth[i] = deptht[i];
    thick[i] = thickt[i];
    outprob1[i] = outprob1t[i];
    outprob2[i] = outprob2t[i];
    type[i] = typet[i];
}

}
