// This function runs the main CPGchron malarkey

#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"use.h"
#include <R_ext/Utils.h>

// Check to see if there are any break requests:
void R_CheckUserInterrupt(void);

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////// Main CPG function //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void cpg(double*thecage,double*thesd,double*thedepths,double*thethick,double*theoutprob1,double*theoutprob2,int*thetypes,char**CALPATH,char**OUTFILE,int*ndets,int*BigCsize,double*LCal,double*HCal,int*m,int*burnin,int*howmany,int*thinby,int*fails)
{

/////////////////////////////////// CONSTANTS ////////////////////////////////////
// Some constants for reading stuff in
// len is calcurve length
int BigCalSize=*BigCsize;
double LowCal=*LCal, HighCal=*HCal;
int IntLowCal = (int)(LowCal*1000);
// Parameters to use on outlier variance
double beta1=2.0,beta2=100.0;
// Let fails be 0 if it succeeds and 1 if it fails
*fails=0;
///////////////////////////// READ IN CALIBRATION CURVE /////////////////////////////

double BigC14[BigCalSize],BigSigma[BigCalSize];

FILE *CalFile;
int i,result;

CalFile = fopen(*CALPATH,"r");

if(CalFile==NULL) {
    error("Error: can't open calibration curve file.\n");
    *fails=1;
    return;
} else {
    Rprintf("Big calibration file opened successfully.\n");

    // Now read in
    for(i=0;i<BigCalSize;i++)
    {
       result=fscanf(CalFile,"%lf",&BigC14[i]);                       
       result=fscanf(CalFile,"%lf",&BigSigma[i]);                                           
    }

    fclose(CalFile);
}	

///////////////////////////// READ IN DETERMINATIONS /////////////////////////////

//enter determinations and their errors - create them as dynamic arrays and enter them
double cage[*ndets],sd[*ndets],depth[*ndets],thick[*ndets],outprob1[*ndets],outprob2[*ndets];
int type[*ndets]; 

Rprintf("========================= \n");
	Rprintf("Data (as read in to 3dp): \n");
    for(i=0;i<*ndets;i++)
    {
       cage[i] = thecage[i];                       
       sd[i] = thesd[i];
       depth[i] = thedepths[i];
       thick[i] = thethick[i];
       outprob1[i] = theoutprob1[i];
       outprob2[i] = theoutprob2[i];
       type[i] = thetypes[i];
	   Rprintf("%3.3lf %3.3lf %3.3lf %3.3lf %3.3lf %3.3lf %i \n",cage[i],sd[i],depth[i],thick[i],outprob1[i],outprob2[i],type[i]);
    }
    Rprintf("End of data. \n");
    Rprintf("========================= \n");


///////////////////// STARTING VALUES ////////////////////////////

// Create starting values for dates
double thetaall[*ndets],shift1[*ndets],shift1new,priorp1[*ndets],shift2[*ndets];
double shift2new,priorp2[*ndets];
double thetanew[*ndets],thetadiff[*ndets-1],thetanewdiff[*ndets-1],badtheta;
double thetanewrat,shift1newrat,psinewrat,meannewrat,shift2newrat;
double currentdepths[*ndets];
int flag1[*ndets],flag1new,flag2[*ndets],flag2new;
int k,badcount;
double p=1.2,mean=5.0,meannew,psi=2.0,psinew; // Starting values
double hi = 10000000; // big number to represent infinity

GetRNGstate();
for(i=0;i<*ndets;i++) 
{
    // Estimate starting values for theta, could do better here with linearinterp
    thetaall[i] = cage[i]+runif(-0.001,0.001);
    flag1[i] = 0;
	shift1[i] = 0;
    flag2[i] = 0;
	shift2[i] = 0;
	priorp1[i] = outprob1[i];
    priorp2[i] = outprob2[i];
}
// depth differences
double depthdiff[*ndets-1];
for(k=0;k<*ndets;k++) currentdepths[k] = depth[k];

int iter;		// iterations loop int;
int q,q1;		// determinations loop int;

// somewhere to store likelihoods
double piytheta[*ndets],pixtheta[*ndets];
double piyflag1[*ndets],pixflag1[*ndets];
double piyshift1[*ndets],pixshift1[*ndets];
double piyflag2[*ndets],pixflag2[*ndets];
double piyshift2[*ndets],pixshift2[*ndets];
double piytwmean=0.0,pixtwmean=0.0;
double piytwpsi,pixtwpsi=0.0;
double baddepth;

// First sort thetas and get differences of all ages
R_qsort(thetaall,1,*ndets);

diff(thetaall,ndets,thetadiff);

///////////////////////////// MCMC PART /////////////////////////////

// Open up the output file
FILE *parameterfile;
parameterfile = fopen(*OUTFILE,"w");

Rprintf("Total number of iterations required: %i \n",*m);
Rprintf("Burn in: %i \n",*burnin);

// Do some timing
clock_t c0, c1, c2;
c0 = clock();

// Let's do some printing to test things that might have gone wrong before the start of iterations
//for(k=0;k<*ndets;k++) Rprintf("k=%i, thetaall = %lf \n",k,thetaall[k]);
//for(k=0;k<*ndets-1;k++) Rprintf("k=%i, thetadiff = %lf \n",k,thetadiff[k]);

// Start iterations here 
for (iter=0;iter<*m;iter++)
{  
	//if(iter==687) {
//		for(k=0;k<*ndets;k++) Rprintf("k=%i, piytheta = %lf \n",k,piytheta[k]);
//		for(k=0;k<*ndets-1;k++) Rprintf("k=%i, pixtheta = %lf \n",k,pixtheta[k]);
//		Rprintf("piytwmean = %lf, pixtwmean = %lf \n",piytwmean,pixtwmean);
//		Rprintf("piytwpsi = %lf, pixtwpsi = %lf \n",piytwpsi,pixtwpsi);
//	} else {
	

    // Get a new seed
    GetRNGstate();

    // Give some updates and some estimated time to finish
    if(iter == *howmany-1) 
    {
        c2 = clock();
        Rprintf("Estimated time to finish is %5.2f minutes or %5.2f seconds \n",(float) (*m/iter)*(c2 - c0)/(60*CLOCKS_PER_SEC),(float) (*m/iter)*(c2 - c0)/(CLOCKS_PER_SEC));
        Rprintf("Iterations so far ... \n");
    }

    // print out some of the values of iter
    if((iter % *howmany == 0) & (iter>=*howmany)) Rprintf("%i \n",iter);

    // Write everything to files
	for(k=0; k<*ndets; k++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", thetaall[k]);
    for(k=0; k<*ndets; k++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", currentdepths[k]);
	for(k=0; k<*ndets; k++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%i ", flag1[k]);
    for(k=0; k<*ndets; k++)	
        if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", shift1[k]);
    for(k=0; k<*ndets; k++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%i ", flag2[k]);
    for(k=0; k<*ndets; k++)	
        if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", shift2[k]);
    if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf %lf \n", mean,psi);

    ///////////////////////////////// DEPTHS ////////////////////////////////////////

    // Use thickness to update the current set of depths used for this particular run
    for(k=0;k<*ndets;k++) {
        currentdepths[k] = runif(depth[k]-thick[k]/2,depth[k]+thick[k]/2);
        // Line to prevent (very rare) duplicates - not used
        //if(k>0) while(currentdepths[k]==currentdepths[k-1]) currentdepths[k] = runif(depth[k]-thick[k]/2,depth[k]+thick[k]/2);
    }
    // And get differences
    diff(currentdepths,ndets,depthdiff);
    // Check for places where the random depths may have swapped everything over
    baddepth = Min(depthdiff,*ndets-1);
    badcount=0;
    while(baddepth<0) {
        // And now re-order all these based on the new currentdepths
        ReOrder(currentdepths,cage,sd,depth,thick,priorp1,priorp2,type,*ndets);
        diff(currentdepths,ndets,depthdiff);
        baddepth = Min(depthdiff,*ndets-1);

        badcount++;
        if(badcount==200) {
            Rprintf("currentdepth cage sd depth thick \n",badcount);
            *fails=1;
            return;
        }
    }

    ///////////////////////////////// THETAS ////////////////////////////////////////

	// Update thetas by looping through each determination
	for(q1=0; q1<*ndets; q1++)
    {       

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

        // All radiocarbon obs updated here
        badtheta = -1.0;
        badcount = 0;
        while(badtheta<0) {
            //sample a new value using a truncated random walk:
            for(i=0;i<*ndets;i++) thetanew[i]= thetaall[i];
    
        	if(q==0)
	       	{
		      	thetanew[0] = truncatedwalk(thetaall[0],0.1*((double)badcount+1),LowCal,thetaall[1]);
			    thetanewrat = truncatedrat(thetaall[0],0.1*((double)badcount+1),LowCal,thetaall[1],thetanew[0]);
		    } else if(q==*ndets-1)
		    {
			    thetanew[*ndets-1] = truncatedwalk(thetaall[*ndets-1],0.1*((double)badcount+1),thetaall[*ndets-2],HighCal);
			    thetanewrat = truncatedrat(thetaall[*ndets-1],0.1*((double)badcount+1),thetaall[*ndets-2],HighCal,thetanew[*ndets-1]);
    		} else {
	       		thetanew[q] = truncatedwalk(thetaall[q],0.1*((double)badcount+1),thetaall[q-1],thetaall[q+1]);		
		      	thetanewrat = truncatedrat(thetaall[q],0.1*((double)badcount+1),thetaall[q-1],thetaall[q+1],thetanew[q]);
		    }
            
            // Difference the new thetas
            diff(thetanew,ndets,thetanewdiff);  
            badtheta = Min(thetadiff,*ndets-1); 
            badcount ++;
            if(badcount==200) {
                Rprintf("Cannot find any satisfactory chronologies. \nCheck the input file for non-consistent depths and ages \n");
                *fails=1;
                return;
            }
        }    
		
	    //calculate old likelihood on first iteration:
        if(iter==0) 
        {
            // Radiocarbon dates 
            if(type[q]==1) {
                pixtheta[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1);
            
                for(i=0;i<*ndets-1;i++) 
                pixtheta[q] += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));
            }

            // Gaussian dates 
            if(type[q]==2) {
                pixtheta[q] = dnorm(thetaall[q],cage[q]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sd[q],1);
            
                for(i=0;i<*ndets-1;i++) 
                pixtheta[q] += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));
            }

            // Uniform dates - tricky this as values outside the uniform range will be 0
            if(type[q]==3) {
                pixtheta[q] = dunif(thetaall[q],cage[q]+flag1[q]*shift1[q]+flag2[q]*shift2[q]-sd[q],cage[q]+flag1[q]*shift1[q]+flag2[q]*shift2[q]+sd[q],1);

                for(i=0;i<*ndets-1;i++) 
                pixtheta[q] += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));
            }

        }

        //calculate new likelihood:
        if(type[q]==1) {
            piytheta[q] = dnorm(cage[q],BigC14[(int)(thetanew[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetanew[q]*1000+0.5)-IntLowCal],2)),1);
            
            for(i=0;i<*ndets-1;i++) {
                piytheta[q] += log(dtweediep1(thetanewdiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1))); 
            }

        } else if(type[q]==2) {
            piytheta[q] = dnorm(thetanew[q],cage[q]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sd[q],1);
            
            for(i=0;i<*ndets-1;i++) 
                piytheta[q] += log(dtweediep1(thetanewdiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));
        } else if(type[q]==3) {
            piytheta[q] = dunif(thetanew[q],cage[q]+flag1[q]*shift1[q]+flag2[q]*shift2[q]-sd[q],cage[q]+flag1[q]*shift1[q]+flag2[q]*shift2[q]+sd[q],1);
            
            for(i=0;i<*ndets-1;i++) 
                piytheta[q] += log(dtweediep1(thetanewdiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));
        }

        //Rprintf("q=%i, pixtheta[q]=%lf, piytheta[q]=%lf \n",q,pixtheta[q],piytheta[q]);
        //if(q==26) Rprintf("th=%lf tp=%i sh1=%lf sh2=%lf \n",thetaall[q],type[q],shift1[q],shift2[q]);
        //if(q==26) Rprintf("fl1=%i \n",flag1[q]);
        //if(q==26) Rprintf("fl2=%i \n",flag2[q]);
        //if(q==26) Rprintf("cage=%lf \n",cage[q]);
        //if(q==26) Rprintf("sd=%lf \n",sd[q]);
        
		//Update the thetas
        thetaall[q] = UpdateMCMC(piytheta[q],pixtheta[q],thetanew[q],thetaall[q],thetanewrat);
	    if(thetaall[q] == thetanew[q]) 
        {
            pixtheta[q] = piytheta[q];
            // Difference the thetas again before re-looping
            diff(thetaall,ndets,thetadiff);
        }

    }	

    // For all dates (radiocarbon or not) the outlier parameters are updated
    for(q1=0; q1<*ndets; q1++)
	{	

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

		// Sample a new flag
		if(flag1[q]==0)
		{
			flag1new = 1;
            if(type[q]==3) flag1new=0;
		} else {
			flag1new = 0;
		}
		
		if(iter==0) 
		{     
		    pixflag1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
                  +dlogbinom(flag1[q],1,priorp1[q]);		
        }
		
        piyflag1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1new*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
              +dlogbinom(flag1new,1,priorp1[q]);		
		
		flag1[q] = (int)UpdateMCMC(piyflag1[q],pixflag1[q],flag1new,flag1[q],1.0);
		if(flag1[q] == flag1new) pixflag1[q] = piyflag1[q];
	}

    for(q1=0;q1<*ndets;q1++)
    {

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

        // Get a new shift
        if(iter==0) 
		{
		       pixshift1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
                    +dnorm(shift1[q],0,sqrt(beta1)*sd[q],1);		
        }
		if(flag1[q]==1)
		{
            shift1new = truncatedwalk(shift1[q],sqrt(beta1)*sd[q],-thetaall[q],26.0-thetaall[q]);
			shift1newrat = truncatedrat(shift1[q],sqrt(beta1)*sd[q],-thetaall[q],26.0-thetaall[q],shift1new);
			
			piyshift1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1new+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
       			  +dnorm(shift1new,0,sqrt(beta1)*sd[q],1);
			
			shift1[q] = UpdateMCMC(piyshift1[q],pixshift1[q],shift1new,shift1[q],shift1newrat);
     		if(shift1[q] == shift1new) pixshift1[q] = piyshift1[q];
		}
        
	}
	//if(iter < 687) {
	// Only going to do this bit for the first 686 iterations - all the rest is done for 687


    for(q1=0; q1<*ndets; q1++)
	{	

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

		// Sample a new big outlier flag
		if(flag2[q]==0)
		{
			flag2new = 1;
            if(type[q]==3) flag2new=0;
		} else {
			flag2new = 0;
		}
		
		if(iter==0) 
		{     
		    pixflag2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
                  +dlogbinom(flag2[q],1,priorp2[q]);		
        }
		
        piyflag2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2new*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
              +dlogbinom(flag2new,1,priorp2[q]);		
		
		flag2[q] = (int)UpdateMCMC(piyflag2[q],pixflag2[q],flag2new,flag2[q],1.0);
		if(flag2[q] == flag2new) pixflag2[q] = piyflag2[q];
	}

    for(q1=0;q1<*ndets;q1++)
    {

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

        // Get a new big outlier shift
        if(iter==0) 
		{
		       pixshift2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
                    +dnorm(shift2[q],0,sqrt(beta2)*sd[q],1);		
        }
		if(flag2[q]==1)
		{
            shift2new = truncatedwalk(shift2[q],sqrt(beta2)*sd[q],-thetaall[q],26.0-thetaall[q]);
			shift2newrat = truncatedrat(shift2[q],sqrt(beta2)*sd[q],-thetaall[q],26.0-thetaall[q],shift2new);
			
			piyshift2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal]+flag1[q]*shift1[q]+flag2[q]*shift2new,sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1)
       			  +dnorm(shift2new,0,sqrt(beta2)*sd[q],1);
			
			shift2[q] = UpdateMCMC(piyshift2[q],pixshift2[q],shift2new,shift2[q],shift2newrat);
     		if(shift2[q] == shift2new) pixshift2[q] = piyshift2[q];
		}
        
	}

    ///////////////////////////////// Tweedie ////////////////////////////////////////

    // Update tweedie parameters
	// First do mean
	meannew = truncatedwalk(mean,0.5,0.0,hi);
	meannewrat = truncatedrat(mean,0.5,0.0,hi,meannew);
	piytwmean = 0.0;
	
	if(iter==0) {
       pixtwmean=0.0;
	   for(i=0;i<*ndets-1;i++) 
          pixtwmean += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));   
       pixtwmean += dgamma(1/mean,0.01,1/0.01,1);
    }
	
	for(i=0;i<*ndets-1;i++) 
		piytwmean += log(dtweediep1(thetadiff[i],p,meannew*depthdiff[i],psi/pow(depthdiff[i],p-1)));
	piytwmean += dgamma(1/meannew,0.01,1/0.01,1);
	
	mean = UpdateMCMC(piytwmean,pixtwmean,meannew,mean,meannewrat);
	if(mean == meannew) pixtwmean = piytwmean;

	// Now update psi
	psinew = truncatedwalk(psi,0.01,0.0,hi);
	psinewrat = truncatedrat(psi,0.01,0.0,hi,psinew);
	piytwpsi = 0.0;

    if(iter==0) {
       pixtwpsi=0.0; 
	   for(i=0;i<*ndets-1;i++) 
          pixtwpsi += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));   
       pixtwpsi += dgamma(1/psi,0.01,1/0.01,1);
    }
    
	for(i=0;i<*ndets-1;i++) 
		piytwpsi += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psinew/pow(depthdiff[i],p-1)));

	piytwpsi += dgamma(1/psinew,0.01,1/0.01,1);
	
	psi = UpdateMCMC(piytwpsi,pixtwpsi,psinew,psi,psinewrat);
	if(psi == psinew) pixtwpsi = piytwpsi;	

    // Sort out the RNG state
    PutRNGstate();

	//}
// End of iterations loop
}

fclose(parameterfile);

c1 = clock();
Rprintf("Completed!\n");
Rprintf("Elapsed time in sec: %5.2f\n",(float) (c1 - c0)/CLOCKS_PER_SEC,2);
Rprintf("Elapsed time in minutes: %5.2f\n",(float) (c1 - c0)/(60ul*CLOCKS_PER_SEC));    
Rprintf("Elapsed time in hours: %5.2f\n",(float) (c1 - c0)/(60ul*60ul*CLOCKS_PER_SEC));

}
