static char mc_version[] = "MATLAB Compiler 1.1 infun";
/*
 *  MATLAB Compiler: 1.1
 *  Date: Jan. 9, 1996
 *  Arguments: -M swvhelp 
 */

 */
#include <math.h>
#include "mex.h"
#include "mcc.h"


void
mexFunction(
    int nlhs_,
    Matrix *plhs_[],
    int nrhs_,
    Matrix *prhs_[]
)
{
   int ci_, i_, j_;
   unsigned flags_;
   Matrix *Mplhs_[32], *Mprhs_[32];
/***************** Compiler Assumptions ****************
 *
 *       I0_         	integer scalar temporary
 *       R0_         	real scalar temporary
 *       R1_         	real scalar temporary
 *       R2_         	real scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       RM1_        	real vector/matrix temporary
 *       counter     	integer scalar
 *       i           	real scalar
 *       j           	real scalar
 *       line        	real vector/matrix
 *       mean        	<function>
 *       mikro       	real scalar
 *       n           	real scalar
 *       nn          	real scalar
 *       points      	real vector/matrix
 *       sd          	real vector/matrix
 *       sum         	<function>
 *       swvhelp     	<function being defined>
 *       ts          	real vector/matrix
 *       tsid        	real vector/matrix
 *       tsmod       	real vector/matrix
 *******************************************************/
   Matrix sd;
   Matrix ts;
   double n = 0.0;
   double i = 0.0;
   int counter = 0;
   double j = 0.0;
   Matrix tsmod;
   double nn = 0.0;
   Matrix points;
   Matrix line;
   double mikro = 0.0;
   Matrix tsid;
   int I0_ = 0;
   Matrix RM0_;
   Matrix RM1_;
   double R0_ = 0.0;
   double R1_ = 0.0;
   double R2_ = 0.0;
   
   mccRealInit(ts);
   mccImport(&ts, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
   n = mccImportReal(0, (nrhs_>1) ? prhs_[1] : 0, "n");
   i = mccImportReal(0, (nrhs_>2) ? prhs_[2] : 0, "i");
   mccRealInit(sd);
   mccRealInit(tsmod);
   mccRealInit(points);
   mccRealInit(line);
   mccRealInit(tsid);
   mccRealInit(RM0_);
   mccRealInit(RM1_);
   
   /* counter=1; */
   counter = 1;
   /* sd=[]; */
   mccCreateEmpty(&sd);
   /* for j=1:2^i:n-i, */
   mccColon(&RM0_, (double)1, pow(2, i), (n - i));
   for (I0_ = 0; I0_ < RM0_.n; ++I0_)
   {
      j = RM0_.pr[I0_];
      /* tsmod=ts(j:j+2^i-1); */
      mccColon2(&RM1_, j, ((j + pow(2, i)) - 1));
      {
         int m_=1, n_=1, cx_ = 0;
         double t_;
         double *p_tsmod;
         int I_tsmod=1;
         double *p_ts;
         double *p_RM1_;
         int I_RM1_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, RM1_.m, RM1_.n, &ts);
         mccAllocateMatrix(&tsmod, m_, n_);
         I_tsmod = (tsmod.m != 1 || tsmod.n != 1);
         p_tsmod = tsmod.pr;
         I_RM1_ = (RM1_.m != 1 || RM1_.n != 1);
         p_RM1_ = RM1_.pr;
         for (j_=0; j_<n_; ++j_)
         {
            for (i_=0; i_<m_; ++i_, p_tsmod+=I_tsmod, p_RM1_+=I_RM1_)
            {
               *p_tsmod = ts.pr[((int)(*p_RM1_ - .5))];
               ;
            }
         }
      }
      tsmod.dmode = mxNUMBER;
      /* %   tsid=bridge(tsmod); */
      
      /* nn=2^i; */
      nn = pow(2, i);
      /* points=(nn-1:-1:0)'; */
      mccColon(&RM1_, (nn - 1), (double) -1, (double)0);
      mccConjTrans(&points, &RM1_);
      /* line=((tsmod(1)-tsmod(nn))*points/(nn-1))+tsmod(nn); */
      R0_ = ((tsmod.pr[(1-1)]) - (tsmod.pr[((int)(nn-.5))]));
      R1_ = (nn - 1);
      R2_ = (tsmod.pr[((int)(nn-.5))]);
      {
         int m_=1, n_=1, cx_ = 0;
         double t_;
         double *p_line;
         int I_line=1;
         double *p_points;
         int I_points=1;
         m_ = mcmCalcResultSize(m_, &n_, points.m, points.n);
         mccAllocateMatrix(&line, m_, n_);
         I_line = (line.m != 1 || line.n != 1);
         p_line = line.pr;
         I_points = (points.m != 1 || points.n != 1);
         p_points = points.pr;
         for (j_=0; j_<n_; ++j_)
         {
            for (i_=0; i_<m_; ++i_, p_line+=I_line, p_points+=I_points)
            {
               *p_line = (((R0_ * *p_points) / (double) R1_) + R2_);
               ;
            }
         }
      }
      line.dmode = mxNUMBER;
      /* tsmod=tsmod-line; */
      {
         int m_=1, n_=1, cx_ = 0;
         double t_;
         double *p_tsmod;
         int I_tsmod=1;
         double *p_1tsmod;
         int I_1tsmod=1;
         double *p_line;
         int I_line=1;
         m_ = mcmCalcResultSize(m_, &n_, tsmod.m, tsmod.n);
         m_ = mcmCalcResultSize(m_, &n_, line.m, line.n);
         mccAllocateMatrix(&tsmod, m_, n_);
         I_tsmod = (tsmod.m != 1 || tsmod.n != 1);
         p_tsmod = tsmod.pr;
         I_1tsmod = (tsmod.m != 1 || tsmod.n != 1);
         p_1tsmod = tsmod.pr;
         I_line = (line.m != 1 || line.n != 1);
         p_line = line.pr;
         for (j_=0; j_<n_; ++j_)
         {
            for (i_=0; i_<m_; ++i_, p_tsmod+=I_tsmod, p_1tsmod+=I_1tsmod, p_line+=I_line)
            {
               *p_tsmod = (*p_1tsmod - *p_line);
               ;
            }
         }
      }
      tsmod.dmode = mxNUMBER;
      
      /* %   sd(counter)=std(tsid); */
      
      /* mikro=mean(tsmod); */
      Mprhs_[0] = &tsmod;
      Mplhs_[0] = 0;
      mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "mean", 15);
      mikro = mccGetReal(0, Mplhs_[ 0 ], " (swvhelp, line 15): mikro");
      /* sd(counter)=(sum((tsmod-mikro).^2)/(nn-1))^0.5; */
      RM1_.dmode = mxNUMBER;
      {
         int m_=1, n_=1, cx_ = 0;
         double t_;
         double *p_RM1_;
         int I_RM1_=1;
         double *p_tsmod;
         int I_tsmod=1;
         m_ = mcmCalcResultSize(m_, &n_, tsmod.m, tsmod.n);
         mccAllocateMatrix(&RM1_, m_, n_);
         I_RM1_ = (RM1_.m != 1 || RM1_.n != 1);
         p_RM1_ = RM1_.pr;
         I_tsmod = (tsmod.m != 1 || tsmod.n != 1);
         p_tsmod = tsmod.pr;
         for (j_=0; j_<n_; ++j_)
         {
            for (i_=0; i_<m_; ++i_, p_RM1_+=I_RM1_, p_tsmod+=I_tsmod)
            {
               *p_RM1_ = mcmRealPowerInt((*p_tsmod - mikro), 2);
               ;
            }
         }
      }
      RM1_.dmode = mxNUMBER;
      sd.pr[(counter-1)] = pow((mccRealVectorSum(&RM1_) / (double) (nn - 1)), 0.5);
      sd.dmode = mxNUMBER;
      
      
      /* tsid=[]; */
      mccCreateEmpty(&tsid);
      /* counter=counter+1; */
      counter = (counter + 1);
      /* end;	 */
   }
   mccReturnFirstValue(&plhs_[0], &sd);
   return;
}
