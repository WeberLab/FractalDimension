function result=acf(data,nlag,nfit,mode)
%===============================================================================
% ACF
%===============================================================================
% 
% Calculates autocorrelation function (ACF) and Hurst coefficient
%
% USAGE:
% result=acf(data,nlag,nfit,mode)
%  
% ARGUMENTS:
% data  : row of data to analyze 
% nlag  : number of lags to calculate, default: 1
% nfit  : number of lags to fit, default: 1. recommended value: 1 to 10.
% mode  : 0: Hurst mode
%         1: autocorrelation function mode
% 
% RETURNED VALUES:
% result: -- Hurst coefficient if mode is 0 or omitted
%         -- autocorrelation function if mode is 1
%
% NEEDS:
% acfopt.m -- rms of difference between ACF and ACF(Hurst) for a given H 
%
% WARNINGS:
% 1, This method is intended only for noise series 
% 2, Fitting a  large number of lags will not give a much better estimate
%    of  the Hurst coefficient than fitting a low  number of lags.  Tests
%    with  synthetic  data  indicate  that  fitting just the first lagged
%    autocorrelation  gives  a marginally lower rms  error for the  Hurst
%    coefficient  for  series of length 64 to 1024 over  a range of Hurst
%    coefficients from 0.1 to 0.9 than fitting the first five lags.
% 3, Due to the limitations of Matlab language in mode 1 the returned ACF
%    has  a  shift  in  its  index.  I.e. result(1)  is  the zeroth  lag,
%    autocorrelation, result(2) is the first lag autocorrelation, etc.
%
% EXAMPLES:
% result=acf(data)               : gives  the Hurst  coefficient of  data
% result=acf(data,nlag,nfit)     : gives  the Hurst  coefficient of  data
%                                  fitted for the first nfit lags 
% result=acf(data,nlag,nfit,1)   : gives the autocorrelation function for
%                                : the first nlag lags
%
%===============================================================================

  
% BACKGROUND:
%   acf  calculates  the autocorrelation function  (ACF) of the  data. An 
%   estimate of the Hurst  coefficient is  also made by fitting  the  ACF
%   with an analytic function  which  depends  on  the Hurst coefficient.
%   In the fitting, each lag is weighted by the inverse of its lag number
%   to  emphasize fitting the lowest lags.  The ACF at j lags is given by
%  
%                SUM from i=1 to N-j (( x(i)-xm)*(x(i+j)-xm) )/(N-j)
%       ACF(j) = ---------------------------------------------------
%                SUM from i=1 to N   (( x(i)-xm)*(x(i  )-xm) )/(N)
%  
%   where j=0 to N-1, N is the number of data values,  x(i) are the  data
%   values,  and xm is the  mean of the x(i)'s.  For a  noise series, the
%   relationship between ACF(j) and the Hurst coefficient, H, is given by
%  
%             ACF(j) = (1/2)( |j+1|^2H -2|j|^2H + |j-1|^2H ).
%  
%   The Hurst estimate is reliable only for noise  series.
%
%   The autocorrelation can be calculated for  as many lags as there  are
%   data  points.  The  analytic  function  which  depends  on the  Hurst
%   coefficient can be required to fit up to the number of lags computed.
%   The  answer  does not vary substantially after the first lag.
%  
%  
%  REFERENCES
%       Bassingthwaighte, J.B. and Beyer, R.P., Fractal  correlation
%       in  heterogeneous  systems (blood flow application). Physica
%       D. vol.53, no.1.  pp. 71-84.  Oct. 1991.
%  
%.......................................................................
% Copyright (C) 1996
% National Simulation Resource, Center for Bioengineering,
% University of Washington, Seattle, WA 98195.
% All rights reserved.
%
% This software may be copied so long as the following copyright notice
% is included:
%
% "This software was developed with support from NIH grant RR-01243"
%
% Please cite this grant in publications for which the software or
% derivatives from it were used and send one reprint to the address
% given above.
%.......................................................................
%  

%===============================================================================
% Starting counter 
%===============================================================================
t0=clock;

%===============================================================================
% Checking arguments   
%===============================================================================
if nargin < 4
   mode = 0;
end

if nargin < 3
   nfit = 1;
end

if nargin < 2
   nlag = 1;
end

n=length(data);
nlag=min(max(1,nlag),n-1);
nfit=min(max(1,nfit),nlag);

%===============================================================================
% Calculating the autocorrelation function   
%===============================================================================
meandata=mean(data);
var=std(data,1)^2;
if var>0
   for tau=0:nlag
      acf(tau+1)=((sum(((data(1:n-tau)-meandata).*(data(1+tau:n)-meandata))))/(n-tau))/var;
   end
else
   acf=zeros(nlag);
end
result=acf;

if mode==1
   elapsed_acf_time=etime(clock,t0)
   return
end

%===============================================================================
% Calculating the Hurst coefficient   
%===============================================================================
clear xhurst options;
xhurst(1)=0.5;%          Starting H value for optimization process
options(1)=0;%           Do not show intermediate steps 
%options(14)=20;%         Maximum number of steps 20 is enough
result=fmins('acfopt',xhurst,options,[],nfit,acf);% optimization  based upon  an 
%                                               analytic relationship of H & acf

%===============================================================================
% Elapsed time
%===============================================================================
elapsed_acf_time=etime(clock,t0)