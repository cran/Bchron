//Calibration script for CPGchron
#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<time.h>

void calibrate(char**CALPATH,char**INFILE,char**OUTFILE,int*ndets)
{


/////////////////////////////////// CONSTANTS ////////////////////////////////////
// Some constants for reading stuff in
// m=no of iterations, len is calcurve length
int m=500000,BigCalSize=26006;
// howmany = how often to print out the number of iterations, thin = how many to thin it by :
int howmany=50000,thinby=5,burnin=50000;

///////////////////////////// READ IN CALIBRATION CURVE /////////////////////////////

double BigC14[BigCalSize],BigSigma[BigCalSize];

FILE *CalFile;
int i;

CalFile = fopen(*CALPATH,"r");

if(CalFile==NULL) {
    Rprintf("Calibration file should be at %s",*CALPATH);
    error("Error: can't open calibration file.");
} else {
    Rprintf("Big calibration file opened successfully.\n");

    // Now read in
    for(i=0;i<BigCalSize;i++)
    {
       fscanf(CalFile,"%lf",&BigC14[i]);                       
       fscanf(CalFile,"%lf",&BigSigma[i]);                                           
    }

    fclose(CalFile);
}	


///////////////////////////// READ IN DETERMINATIONS /////////////////////////////

//enter determinations and their errors - create them as dynamic arrays and enter them
//from a separate file using same method as cal curve:
char labcode[*ndets][50];
double cage[*ndets],sd[*ndets],depth[*ndets],thick[*ndets],outprob1[*ndets],outprob2[*ndets];
int type[*ndets];
 
FILE *dets;

double numb1[*ndets],numb2[*ndets],numb3[*ndets],numb4[*ndets],numb5[*ndets],numb6[*ndets];
int numb7[*ndets];

dets = fopen(*INFILE,"r");

if(dets==NULL) {
    error("Error: can't open determinations file.");
} else {
    Rprintf("Determinations file opened successfully.\n");

    // First get rid of header
    char temp[100];
    fgets(temp,100,dets);

    // Now read in as ints and then loop again to convert to double
    for(i=0;i<*ndets;i++)
    {
       fscanf(dets,"%s",labcode[i]);   
       fscanf(dets,"%lf",&numb1[i]);                       
       fscanf(dets,"%lf",&numb2[i]);                       
       fscanf(dets,"%lf",&numb3[i]);                       
       fscanf(dets,"%lf",&numb4[i]);                       
       fscanf(dets,"%lf",&numb5[i]);
       fscanf(dets,"%lf",&numb6[i]); 
       fscanf(dets,"%i",&numb7[i]); 
    }

    for(i=0;i<*ndets;i++)
    {
       cage[i] = (double)numb1[i]/1000;                       
       sd[i] = (double)numb2[i]/1000;
       depth[i] = (double)numb3[i]/100;
       thick[i] = (double)numb4[i]/100;
       outprob1[i] = (double)numb5[i];
       outprob2[i] = (double)numb6[i];
       type[i] = numb7[i];
    }
    
    Rprintf("Determinations read successfully.\n");
    
    fclose(dets);
}



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
		for(i=0;i<*ndets;i++) thetanew[i]= thetaall[i];
		thetanew[q] = rnorm(thetaall[q],0.1);
        
        // Stop it from choosing bad values outside the range of BigCal
        while((thetanew[q]< -0.005) | (thetanew[q] > 26)) thetanew[q] = rnorm(thetaall[q],0.1);

        //if(iter % thinby == 0 &&  iter > burnin) Rprintf("%lf \n",thetanew[2]);     

		//calculate old likelihood on first iteration:
        if(iter==0) 
        {
           pixtheta[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1);
        }

        //calculate new likelihood:
        piytheta[q] = dnorm(cage[q],BigC14[(int)(thetanew[q]*1000+0.5)+5],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetanew[q]*1000+0.5)+5],2)),1);
     
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
