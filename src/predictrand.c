//This file runs a random prediction stage of the CPGchron method

#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<time.h>
#include"use.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////// Random Predict function ////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void predictrand(char**PARFILE,char**DETSFILE,char**OUTFILE,double *lowddepths,double *highddepths,int *nddepthints, int*ndets,int*numchrons,double*Present,char**OUTLIERFILE)
{

///////////////////////////// READ IN DETERMINATIONS /////////////////////////////

//enter determinations and their errors - create them as dynamic arrays and enter them
//from a separate file using same method as cal curve:
char labcode[*ndets][50];
double cage[*ndets],sd[*ndets],depth[*ndets],thick[*ndets],outprob1[*ndets],outprob2[*ndets];
int type[*ndets]; 
double currentddepth[1];

FILE *dets;

double numb1[*ndets],numb2[*ndets],numb3[*ndets],numb4[*ndets],numb5[*ndets],numb6[*ndets];
int numb7[*ndets];
int i;

dets = fopen(*DETSFILE,"r");

if(dets==NULL) {
    error("Error: can't open determinations file.\n");
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
       type[i] = (int)numb7[i];
    }
    
    Rprintf("Determinations read successfully.\n");
    
    fclose(dets);
}

///////////////////// STARTING VALUES ////////////////////////////

// Set up the value of p
double p = 1.2;

//////////////////////// READ IN PARAMETER AND START PREDICTION //////////////////////////

// Create arrays to store everything
double thetas[*ndets],shift1[*ndets],shift2[*ndets],mydepths[*ndets];
double mean,psi;
int flag1[*ndets],flag2[*ndets];
double PredEst,alphaT,lambdaT,betaT;
int OutlierSum1[*ndets],OutlierSum2[*ndets];
int Nd;
int j,k,K=1000,wrong=0,count=0;
int stopper, counter;
double unitemp;

// Set PredEst to zero
PredEst = 0.0;

// Set the OutlierSums to zero
for(k=0;k<*ndets;k++) OutlierSum1[k] = 0.0;
for(k=0;k<*ndets;k++) OutlierSum2[k] = 0.0;

// This is a tricky one - need to get the design depths in each segment
int lencurrentdepths;

// Calculate the differences of the depths
double depthdiff[*ndets-1];

FILE *pars,*chrons;

// Cleverly, I'm going to read it in one line at a time so that I don't have to store too much
pars = fopen(*PARFILE,"r");
chrons = fopen(*OUTFILE,"w");

if(pars==NULL) {
    error("Error: can't open file of parameters.\n");
} else {
    Rprintf("Parameters file opened successfully.\n");

    // Read in pars one line at a time
    for(i=0;i<*numchrons;i++)
    {

       for(j=0;j<*ndets;j++)
       {  
         fscanf(pars,"%lf",&thetas[j]);                       
       } 
       for(j=0;j<*ndets;j++)
       {  
         fscanf(pars,"%lf",&mydepths[j]);                       
       } 
       for(j=0;j<*ndets;j++)
       {  
         fscanf(pars,"%i",&flag1[j]);                       
       } 
       for(j=0;j<*ndets;j++)
       {  
         fscanf(pars,"%lf",&shift1[j]);                       
       } 
       for(j=0;j<*ndets;j++)
       {  
         fscanf(pars,"%i",&flag2[j]);                       
       } 
       for(j=0;j<*ndets;j++)
       {  
         fscanf(pars,"%lf",&shift2[j]);                       
       } 
       fscanf(pars,"%lf",&mean);                       
       fscanf(pars,"%lf",&psi);                       

    alphaT = (2-p)/(p-1);
    lambdaT = pow(mean,(2-p))/(psi*(2-p));
    betaT = 1/(psi*(p-1)*pow(mean,(p-1)));
    diff(mydepths,ndets,depthdiff);

    
    for(k=0;k<*ndets;k++)
    {
        OutlierSum1[k] +=flag1[k];
        OutlierSum2[k] +=flag2[k];
    }

    if(i % 100==0) Rprintf("%i \n", *numchrons-i);

    // Get a new seed
    GetRNGstate();

    // Get the currentddepth
    unitemp = runif(0,1);
    stopper = -1;
    counter = 0;
    while(stopper<0) {
        if(unitemp<(double)(counter+1)/(*nddepthints)) {
            currentddepth[0] = runif(lowddepths[counter],highddepths[counter]);
            stopper = 1;
        }
        counter++;
    }

    // Now loop through depth to get a set of values for each of the interpolates
    for(j=0;j<*ndets-1;j++) 
    {

        Nd = rpois(lambdaT*depthdiff[j]);  

        // Get the current depths and find the indices at which they are at          
        lencurrentdepths = GetLengthCurrentDepths(mydepths[j],mydepths[j+1],currentddepth,1);

        if(lencurrentdepths>0) {

            // Generate a load of runif/rexps which give the depth cutoffs for each section
            double Tempexp[Nd];
            int Tempexp2[Nd];
            if(Nd>0) 
            {
            for(k=0;k<Nd;k++)
                Tempexp[k] = runif(mydepths[j],mydepths[j+1]);
            
            for(k=0;k<Nd;k++) Tempexp2[k] = (int)(Tempexp[k]*100000);
            qsort(Tempexp2,Nd,sizeof(int),compare);
            for(k=0;k<Nd;k++) Tempexp[k] = (double)(Tempexp2[k])/100000;
    
            }
           
            // Now generate the gamma bits
            double Temp[Nd+1],Temp2,Temp3[Nd+1],Temp4[Nd+1];
            Temp2 = 0.0;
            for(k=0;k<Nd+1;k++)
            {
                Temp[k] = rgamma(alphaT,1/betaT);
                Temp2 += Temp[k];
            }
            for(k=0;k<Nd+1;k++) Temp3[k] = (thetas[j+1]-thetas[j])*Temp[k]/Temp2;
            Temp4[0] = Temp3[0];
            for(k=1;k<Nd+1;k++) Temp4[k] = Temp4[k-1]+Temp3[k];
    
            // Now do some linear interpolation - setup
            double xinterp[Nd+2],yinterp[Nd+2];
            xinterp[0] = mydepths[j];
            xinterp[Nd+1] = mydepths[j+1];
            if(Nd>0) for(k=1;k<Nd+1;k++) xinterp[k] = Tempexp[k-1];
            yinterp[0] = thetas[j];
            for(k=1;k<Nd+2;k++) yinterp[k] = thetas[j]+Temp4[k-1];
    
            // Interpolation steps
            PredEst = linearinterp(Nd+2, currentddepth[0], xinterp, yinterp);
     
        }
    }

    // Now need to extrapolate beyond the edge of depth
    
    // First upwards
    // Find the depths that are above the first depth
    lencurrentdepths = GetLengthCurrentDepths(0,mydepths[0],currentddepth,1);

    if(lencurrentdepths>0)
    {
        // Loop so that none of the PredDates are in the future.
        wrong = 1;
        count = 0;
        double Temp,Tempexp;
        double Temp2[K],Tempexp2[K],xinterp[K],yinterp[K];
        while(wrong>0)
        {
            count++;
            for(k=0;k<K;k++) 
            {
                Tempexp = rexp(1/lambdaT);
                if(k==0) {
                    Tempexp2[0] = Tempexp;
                } else {
                    Tempexp2[k] = Tempexp2[k-1]+Tempexp;
                }
            }
            for(k=0;k<K;k++) Tempexp2[k] = mydepths[0]-Tempexp2[k];
        
            for(k=0;k<K;k++) 
            {
                Temp = rgamma(alphaT,1/betaT);
                if(k==0) {
                    Temp2[0] = Temp;
                } else {
                    Temp2[k] = Temp2[k-1]+Temp;
                }
            }
            for(k=0;k<K;k++) Temp2[k] = thetas[0]-Temp2[k];
        
            // Interpolation steps
            xinterp[0] = mydepths[0];
            for(k=1;k<K;k++) xinterp[k] = Tempexp2[k];
            yinterp[0] = thetas[0];
            for(k=1;k<K;k++) yinterp[k] = Temp2[k];
    
            PredEst = linearinterp(K, currentddepth[0], xinterp, yinterp);
    
            if(PredEst<*Present) wrong = -1;
        }
            

    }
    
    // Now interpolate below
    // Find the depths that are above the first depth
    lencurrentdepths = GetLengthCurrentDepths(mydepths[*ndets-1],currentddepth[0],currentddepth,1);
    if(lencurrentdepths>0)
    {
        //  Now do the creating again
        double Tempexp,Tempexp2[K];
        
        Tempexp2[0] = 0.0;
        for(k=1;k<K;k++) 
        {
            Tempexp = rexp(1/lambdaT);
            Tempexp2[k] = Tempexp2[k-1]+Tempexp;
        }
        for(k=0;k<K;k++) Tempexp2[k] = mydepths[*ndets-1]+Tempexp2[k];
        
        double Temp,Temp2[K];
        Temp2[0] = 0.0;
        for(k=1;k<K;k++) 
        {
            Temp = rgamma(alphaT,1/betaT);
            Temp2[k] = Temp2[k-1]+Temp;
        }
        for(k=0;k<K;k++) Temp2[k] = thetas[*ndets-1]+Temp2[k];
        
        // Interpolation steps
        double xinterp[K],yinterp[K];
        for(k=0;k<K;k++) xinterp[k] = Tempexp2[k];
        for(k=0;k<K;k++) yinterp[k] = Temp2[k];
    
        PredEst = linearinterp(K, currentddepth[0], xinterp, yinterp);
    
    }   
 
    // Get a new seed
    PutRNGstate();

    /////////////////////////////////////////////////////////////////////////////
   
    
    // Now write out to file
	fprintf(chrons,"%lf \n", PredEst);
    

    }

    fclose(pars);
    fclose(chrons);

    // Write out the outliers to a file
    FILE *outliers;
    outliers = fopen(*OUTLIERFILE,"w");
    fprintf(outliers,"Labcode Prob1 Prob2 \n");
    for(k=0; k<*ndets; k++)	
    {
        fprintf(outliers,"%s ",labcode[k]);
		fprintf(outliers,"%5.3lf %5.3lf \n", (double)OutlierSum1[k]/ (double)*numchrons,(double)OutlierSum2[k]/ (double)*numchrons);
    }
    fclose(outliers);



}

}
