#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <alloc.h>

/* calculate V(j) complex matrix with dhmhelp */
void dhmhelp(double *Vr, double *Vi, double *s, double *w, long n, long M)
{ long i, j;
double sqroh, temp;

sqroh=pow(0.5,0.5);
for (j=0; j<M; j++) {
	temp=*(s+j);
	*(s+j)=pow(temp,0.5);}
*Vr=*s * *w;
*(Vr+n)=*(s+n) * *(w+M-1);
*Vi=0;
*(Vi+n)=0;

i=1;
for (j=1; j<n; j++) {
	*(Vr+j)=sqroh * *(s+j) * *(w+i);
   i++;
   *(Vi+j)=*(w+i);
   i++;}
return;
}

/* The gateway function */
void mexFunction (
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
   {long n, M, rows;
   double *Vr, *Vi, *s, *w;

   /* Check for proper number of arguments. */
   if (nrhs !=3) {
   	mexErrMsgTxt ("Three inputs required.");}
   else if (nlhs !=1) {
      mexErrMsgTxt ("Only one output is allowed.");}

   /* Check to make sure the third input arguments is a scalar. */
   if (mxIsComplex(prhs[2]) || mxGetN(prhs[2]) *  mxGetM(prhs[2]) != 1) {
   	mexErrMsgTxt ("Input n must be scalar.");}

   /* Get the scalar inputs. */
   n=mxGetScalar(prhs[2]);

   /* Create a pointer to a copy of the first and second input matrix ts. */
   s=mxGetPr(prhs[0]);
   w=mxGetPr(prhs[1]);

   /* Assign the pointer to the output variable. */
   rows=n+1;
   plhs[0]=mxCreateDoubleMatrix(rows,1,mxCOMPLEX);
   Vr=mxGetPr(plhs[0]);
   Vi=mxGetPi(plhs[0]);

   /* Call the dhmhelp c subroutine. */
   M=mxGetN(prhs[1])*mxGetM(prhs[1]);
 	dhmhelp(Vr, Vi, s, w, n, M);}
