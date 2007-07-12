//This file runs the Prediction stage of the CPGchron method

// This function runs the main CPGchron malarkey

#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<time.h>
#include"use.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////// Main Predict function //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void predict(char**PARFILE,char**DETSFILE,char**OUTFILE,int*ndets,char**DDEPTHFILE,int*lenddepths,int*numchrons,double*Present,char**OUTLIERFILE)
{

///////////////////////////// READ IN DETERMINATIONS /////////////////////////////

//enter determinations and their errors - create them as dynamic arrays and enter them
//from a separate file using same method as cal curve:
char labcode[*ndets][50];
double cage[*ndets],sd[*ndets],depth[*ndets],thick[*ndets],outprob1[*ndets],outprob2[*ndets];
int type[*ndets]; 

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

//////////////////////// READ IN DDEPTHS //////////////////////////

// Read in the 14C depths and the design depths
double ddepths[*lenddepths];

FILE *ddfile;

// Open design depths file first

ddfile = fopen(*DDEPTHFILE,"r");

if(ddfile==NULL) {
    error("Error: can't open design depths file.\n");
} else {
    Rprintf("Design depths file opened successfully.\n");

    // Now read in 
    for(i=0;i<*lenddepths;i++)
    {
       fscanf(ddfile,"%lf",&ddepths[i]);  
       ddepths[i] = ddepths[i]/100;  
    }

    Rprintf("Design depths read successfully.\n");

    fclose(ddfile);
}

///////////////////// STARTING VALUES ////////////////////////////

// Set up the value of p
double p = 1.2;

//////////////////////// READ IN PARAMETER AND START PREDICTION //////////////////////////

// Create arrays to store everything
double thetas[*ndets],shift1[*ndets],shift2[*ndets],mydepths[*ndets];
double mean,psi;
int flag1[*ndets],flag2[*ndets];
double PredEst[*lenddepths],alphaT,lambdaT,betaT;
int OutlierSum1[*ndets],OutlierSum2[*ndets];
int Nd;
int j,k,m,K=1000,wrong=0,count=0;

// Set PredEsts to zero
for(k=0;k<*lenddepths;k++) PredEst[k] = 0.0;

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

    // Now loop through depth to get a set of values for each of the interpolates
    for(j=0;j<*ndets-1;j++) 
    {

        Nd = rpois(lambdaT*depthdiff[j]);  

        // Get the current depths and find the indices at which they are at          
        lencurrentdepths = GetLengthCurrentDepths(mydepths[j],mydepths[j+1],ddepths,*lenddepths);
        double currentdepths[lencurrentdepths];
        int currentdepthrows[lencurrentdepths];
        GetCurrentDepths(mydepths[j],mydepths[j+1],ddepths,*lenddepths,currentdepths);
        GetCurrentDepthRows(mydepths[j],mydepths[j+1],ddepths,*lenddepths,currentdepthrows);     

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
        double PredDates[lencurrentdepths];
        double xinterp[Nd+2],yinterp[Nd+2];
        xinterp[0] = mydepths[j];
        xinterp[Nd+1] = mydepths[j+1];
        if(Nd>0) for(k=1;k<Nd+1;k++) xinterp[k] = Tempexp[k-1];
        yinterp[0] = thetas[j];
        for(k=1;k<Nd+2;k++) yinterp[k] = thetas[j]+Temp4[k-1];

        // Interpolation steps
        for(k=0;k<lencurrentdepths;k++) PredDates[k] = linearinterp(Nd+2, currentdepths[k], xinterp, yinterp);
 
        // Now fill in the blank spaces in the bigger PredEst
        for(k=0;k<*lenddepths;k++) 
        {
            for(m=0;m<lencurrentdepths;m++)
            {
                if(k==currentdepthrows[m]) PredEst[currentdepthrows[m]] = PredDates[m];
            }
        }     

    }

    // Now need to extrapolate beyond the edge of depth
    
    // First upwards
    // Find the depths that are above the first depth
    lencurrentdepths = GetLengthCurrentDepths(0,mydepths[0],ddepths,*lenddepths);

    if(lencurrentdepths>0)
    {
    double currentdepths[lencurrentdepths];
    int currentdepthrows[lencurrentdepths];
    GetCurrentDepths(0,mydepths[0],ddepths,*lenddepths,currentdepths);
    GetCurrentDepthRows(0,mydepths[0],ddepths,*lenddepths,currentdepthrows); 

    //  Now do the creating again
    double Tempexp=0.0,Tempexp2[K];
    
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
   
    double Temp,Temp2[K];
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
    double PredDates[lencurrentdepths];
    double xinterp[K],yinterp[K];
    xinterp[0] = mydepths[0];
    for(k=1;k<K;k++) xinterp[k] = Tempexp2[k];
    yinterp[0] = thetas[0];
    for(k=1;k<K;k++) yinterp[k] = Temp2[k];

    for(k=0;k<lencurrentdepths;k++) PredDates[k] = linearinterp(K, currentdepths[k], xinterp, yinterp);

    // Re-loop if any of the PredDates are in the future.
    wrong = 0;
    for(k=0;k<lencurrentdepths;k++) if(PredDates[k]<*Present) wrong += 1;

    while(wrong>0)
    {
        count += 1;
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

        for(k=0;k<lencurrentdepths;k++) PredDates[k] = linearinterp(K, currentdepths[k], xinterp, yinterp);

        wrong = 0;
        for(k=0;k<lencurrentdepths;k++) if(PredDates[k]<*Present)
        {
            wrong += 1;
        }
        if(wrong == 0) 
        { 
            count = 0;       
        }
 
        if(count==500) 
        {
            for(k=0;k<lencurrentdepths;k++) if(PredDates[k]<*Present) PredDates[k] = *Present;
            wrong = 0;
            count = 0;
        }
       
    }
  
    // and write to the bigger array
    for(k=0;k<*lenddepths;k++) 
    {
        for(m=0;m<lencurrentdepths;m++)
        {
            if(k==currentdepthrows[m]) PredEst[currentdepthrows[m]] = PredDates[m];
        }
    }   
    }
    
    // Now interpolate below
    // Find the depths that are above the first depth
    lencurrentdepths = GetLengthCurrentDepths(mydepths[*ndets-1],ddepths[*lenddepths-1],ddepths,*lenddepths);
    if(lencurrentdepths>0)
    {
    double currentdepths[lencurrentdepths];
    int currentdepthrows[lencurrentdepths];
    GetCurrentDepths(mydepths[*ndets-1],ddepths[*lenddepths-1],ddepths,*lenddepths,currentdepths);
    GetCurrentDepthRows(mydepths[*ndets-1],ddepths[*lenddepths-1],ddepths,*lenddepths,currentdepthrows); 

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
    double PredDates[lencurrentdepths];
    double xinterp[K],yinterp[K];
    //xinterp[0] = depth[ndets];
    for(k=0;k<K;k++) xinterp[k] = Tempexp2[k];
    //yinterp[0] = thetas[ndets];
    for(k=0;k<K;k++) yinterp[k] = Temp2[k];

    for(k=0;k<lencurrentdepths;k++) PredDates[k] = linearinterp(K, currentdepths[k], xinterp, yinterp);

    for(k=0;k<*lenddepths;k++) 
    {
        for(m=0;m<lencurrentdepths;m++)
        {
            if(k==currentdepthrows[m]) PredEst[currentdepthrows[m]] = PredDates[m];
        }
    }   
    
    }

    // Get a new seed
    PutRNGstate();

    /////////////////////////////////////////////////////////////////////////////
   
    
    // Now write out to file
    for(k=0; k<*lenddepths; k++)	
    {
		fprintf(chrons,"%lf ", PredEst[k]);
    }
    fprintf(chrons,"\n");
    

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
