function result = dfa(data,slowlen,midlen,fastlen,step)
%===============================================================================
% DFA
%===============================================================================
%
% Compute DFA in the manner described by C-K Peng,
% in C-K Peng, S Havlin, HE Stanley, AL Goldberger (1995)
% "Detrended Fluctuation Analysis" Chaos 5(1):82-87
%
% USAGE:
% result = dfa(vals,slowlen,midlen,fastlen,step)
%
% ARGUMENTS:
% data     : the data
% slowlen  : a slow time scale, 64 by default
% midlen   : a medium time scale, 16 by default
% fastlen  : a fast time scale, 4 by default
% step     : step between the size of following boxes
%
% RETURNED VALUES:
% result(1): Long range scaling exponent (alpha) 
% result(2): Long range spectral index (beta) 
% result(3): Long range Hurst coefficient (H) for fractional Gaussian noise 
%            (fGn) series
% result(4): Long range Hurst coefficient (H) for fractional Brownian motion
%            (fBm) series
% result(5): Short range scaling exponent (alpha) 
% result(6): Short range spectral index (beta) 
% result(7): Short range Hurst coefficient (H) for fractional Gaussian noise
%            (fGn) series
% result(8): Short range Hurst coefficient (H) for fractional Brownian motion
%            (fBm) series
%
% WRITTEN BY:
% Lajos R. Kozak
%
% Program based on D. Kaplan's Matlab version of the DFA algorithm
% 
% Modifications:
% -- faster data integration
% -- spectral index & Hurst exponent as returned values
% -- user definable stepsize
% -- different parts integrated into one program
%
% NEEDS:
% detrend.m -- from the Signal Processing Toolbox
%
%===============================================================================

%===============================================================================
% Starting counter 
%===============================================================================
t0=clock;

%===============================================================================
% Checking arguments   
%===============================================================================
if nargin < 5
   step = sqrt(2);
end

if nargin < 4
   fastlen = 4;
end

if nargin < 3
   midlen = 16;
end

if nargin < 2
   slowlen = 64;
end

%===============================================================================
% Integrating data
%===============================================================================
intdata=cumsum(data-mean(data));

num=0;
resid=0;
N = length(intdata);
i = fastlen;
j = 1;

%===============================================================================
% Doing detrending and calculating RMS for different boxsizes
%===============================================================================
while floor(i) <= slowlen
   
   boxlength=floor(i);
   sumsq = 0;
   datausedcount=0;
   for k = 0:boxlength:N-1
      if (k+boxlength > N)
         break
      end
      datausedcount = k+boxlength;
      ypoint = intdata(k+1:k+boxlength);
      %=========================================================================
      % Detrending
      %=========================================================================
      ypoint = detrend(ypoint); 
      sumsq = sumsq+sum(ypoint.^2);
   end
   
   %============================================================================
   % Calculating RMS
   %============================================================================
   resid(j) = sqrt(sumsq/datausedcount);
   num(j) = (floor(i));
   i = i*step;
   j = j+1;

end

%===============================================================================
% Line fitting (long range)
%===============================================================================
goods = find( num >= midlen & num <= slowlen);
pl = polyfit( log(num(goods)), log(resid(goods)), 1 );

%===============================================================================
% Line fitting (short range)
%===============================================================================
goods = find( num >= fastlen & num <= midlen);
ps = polyfit( log(num(goods)), log(resid(goods)), 1 );

%===============================================================================
% Calculating results
%===============================================================================
result(1) = pl(1);%            alpha         
result(2) = 2*result(1)-1;%    beta
result(3) = (result(2)+1)/2;%  Hurst (fGn)
result(4) = (result(2)-1)/2;%  Hurst (fBm)
result(5) = ps(1);%            alpha         
result(6) = 2*result(5)-1;%    beta
result(7) = (result(6)+1)/2;%  Hurst (fGn)
result(8) = (result(6)-1)/2;%  Hurst (fBm)

%===============================================================================
% Elapsed time
%===============================================================================
elapsed_dfa_time=etime(clock,t0)