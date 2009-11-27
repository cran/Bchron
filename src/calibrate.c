//Calibration script for CPGchron
#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<time.h>

void calibrate(double*c14dates,double*c14errors,int*datetypes,char**CALPATH,char**OUTFILE,int*ndets,int*BigCSize,double*LCal,double*HCal,int*iterat,int*burn,int*thin,int*howmny)
{


/////////////////////////////////// CONSTANTS ////////////////////////////////////
// Some constants for reading stuff in
// m=no of iterations, len is calcurve length
int m=*iterat,BigCalSize=*BigCSize;
double LowCal=*LCal, HighCal=*HCal;
int IntLowCal = (int)(LowCal*1000);
// howmany = how often to print out the number of iterations, thin = how many to thin it by :
int howmany=*howmny,thinby=*thin,burnin=*burn;

///////////////////////////// READ IN CALIBRATION CURVE /////////////////////////////

double BigC14[BigCalSize],BigSigma[BigCalSize];

FILE *CalFile;
int i;

CalFile = fopen(*CALPATH,"r");

if(CalFile==NULL) {
    Rprintf("Calibration file should be at %s",*CALPATH);
    error("Error: can't open calibration file.");
} else {
    Rprintf("Calibration file opened successfully.\n");

    // Now read in
    for(i=0;i<BigCalSize;i++)
    {
       fscanf(CalFile,"%lf",&BigC14[i]);                       
       fscanf(CalFile,"%lf",&BigSigma[i]);                                           
    }

    fclose(CalFile);
}	


///////////////////////////// READ IN DETERMINATIONS /////////////////////////////

double cage[*ndets],sd[*ndets];
int type[*ndets];

Rprintf("================ \n");
Rprintf("Date  Error Type \n");
for(i=0;i<*ndets;i++) {
    cage[i] = c14dates[i];                       
	sd[i] = c14errors[i];
	type[i] = datetypes[i];
    Rprintf("%4.3lf %4.3lf %i \n",cage[i],sd[i],type[i]);
}
Rprintf("================ \n");

///////////////////////////// CALIBRATION MCMC //////////////////////////    

// Do some timing
clock_t c0, c1, c2;
c0 = clock();

// Declare variables
double thetaall[*ndets],thetanew[*ndets];
double piytheta[*ndets],pixtheta[*ndets];
int iter,q;
double U;

// Setup initial thetas - some more clever things could be done here
for(i=0;i<*ndets;i++) thetaall[i] = cage[i];

// Open up the output file
FILE *parameterfile;
parameterfile = fopen(*OUTFILE,"w");

Rprintf("Total number of iterations required: %i \n",m);
Rprintf("Burn-in size: %i \n",burnin);
Rprintf("Thinning by: %i \n",thinby);

// Start iterations here 
for (iter=0;iter<m;iter++)
{

    // print out some of the values of iter
    if(iter % howmany == 0 && iter>1)	Rprintf("%i \n",iter);
    
    // Give some update and some estimated time to finish
    if(iter == 10000) 
    {
        c2 = clock();
        Rprintf("Estimated time to finish is %5.2f minutes or %5.2f seconds \n",(float) (m/iter)*(c2 - c0)/(60*CLOCKS_PER_SEC),(float) (m/iter)*(c2 - c0)/(CLOCKS_PER_SEC));
        Rprintf("Iterations so far ... \n");
    }

	// Write everything to files
    if(iter % thinby == 0 && iter > burnin) 
    {
    	for(q=0; q<*ndets; q++)	fprintf(parameterfile,"%lf ", thetaall[q]);
        fprintf(parameterfile,"\n");
    }
 
    // Update thetas by looping through each determination
	for(q=0; q<*ndets; q++)
    { 
        
        // Get a new seed
        GetRNGstate(); 

        // Only update radiocarbon dates that are of type 1, otherwise use raw date.
        if(type[q]==2) {
            if(sd[q]==0) {
                thetaall[q] = cage[q];
            } else {
                thetaall[q] = rnorm(cage[q],sd[q]);
            }
        } 
        else if(type[q]==3) {
            if(sd[q]==0) {
                thetaall[q] = cage[q];
            } else {
                thetaall[q] = runif(cage[q]-sd[q],cage[q]+sd[q]);
            }
        } else {

     	//sample a new value from a distribution using function from random.c:
		thetanew[q] = rnorm(thetaall[q],runif(0.0,1.0));
        
        // Stop it from choosing bad values outside the range of BigCal
        if(type[q]==1) while((thetanew[q]< LowCal) | (thetanew[q] > HighCal)) thetanew[q] = rnorm(thetaall[q],0.1);

        //if(iter % thinby == 0 &&  iter > burnin) Rprintf("%lf \n",thetanew[2]);     

		//calculate old likelihood on first iteration:
        if(iter==0) 
        {
           pixtheta[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)-IntLowCal],2)),1);
        }

        //calculate new likelihood:
        piytheta[q] = dnorm(cage[q],BigC14[(int)(thetanew[q]*1000+0.5)-IntLowCal],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetanew[q]*1000+0.5)-IntLowCal],2)),1);
        /*if(q==18) {
            Rprintf("cage[q]=%lf, sd[q]=%lf, LowCal=%lf, HighCal=%lf, IntLowCal=%i \n",cage[q],sd[q],LowCal,HighCal,IntLowCal);
            Rprintf("thetaall[q]=%lf, thetanew[q]=%lf,BigC14[thetaall[q]]=%lf,BigC14[thetanew[q]]=%lf \n",thetaall[q],thetanew[q],BigC14[(int)(thetaall[q]*1000+0.5)-IntLowCal],BigC14[(int)(thetanew[q]*1000+0.5)-IntLowCal]);
            Rprintf("pixtheta[q]=%lf, piytheta[q]=%lf \n",pixtheta[q],piytheta[q]);
        }*/
       

		//Update the thetas
        U = runif(0.0,1.0);
        if(U<exp(piytheta[q]-pixtheta[q])) thetaall[q] = thetanew[q];
		if(thetaall[q] == thetanew[q]) pixtheta[q] = piytheta[q];

        }

        PutRNGstate();
   
	}     

} 
  
fclose(parameterfile);

c1 = clock();
Rprintf("Completed!\n");
Rprintf("Elapsed time in sec: %5.2f\n",(float) (c1 - c0)/CLOCKS_PER_SEC,2);
Rprintf("Elapsed time in minutes: %5.2f\n",(float) (c1 - c0)/(60*CLOCKS_PER_SEC));    
Rprintf("Elapsed time in hours: %5.2f\n",(float) (c1 - c0)/(60*60*CLOCKS_PER_SEC));


}