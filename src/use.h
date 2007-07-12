// Header file for use.c

int diff(double *arr,int *len,double *retarr);
int compare (const void * a, const void * b);
int GetLengthCurrentDepths(double depthlow,double depthhigh,double *ddepth,int length);
int GetCurrentDepths(double depthlow,double depthhigh,double *ddepth,int length,double *current);
int GetCurrentDepthRows(double depthlow,double depthhigh,double *ddepth,int length,int *rows);
double linearinterp(int n, double newx, double *a, double *b);
double Max2 (double a, double b);
double rtruncn (double a, double b);
double truncatedwalk (double old, double sd, double low, double high);
double truncatedrat (double old, double sd, double low, double high, double newvalue);
double Max(double *Numbers, int Count);
double Min(double *Numbers, int Count);
int seq(double from,double to,double len,double *sequence);
double dtweedielogwsmallp(double y, double phi, double power);
double dtweedieseriessmallp(double power,double y, double mu, double phi);
double dtweediep1(double y, double power, double mu, double phi);
double UpdateMCMC(double newloglik,double oldloglik,double newval,double oldval,double rat);
int fact(int number);
double dlogbinom(int x,int n, double p);
int compare_doubles (const double *a, const double *b);
void ReOrder(double *currentdepths,double *cage,double *sd,double *depth,double *thick, double *outprob1, double *outprob2, int *type, int ndets);
