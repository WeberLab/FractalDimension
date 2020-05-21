function rms=acfopt(x,nfit,auto);
%===============================================================================
% ACFOPT
%===============================================================================
% Calculate rms of weighted differenced between the ideal ACF of a known Hurst
% coefficient and the given ACF in order to determine H from a known ACF 
%
% USAGE:
% rms=acfopt(x,nfit,auto)
%  
% ARGUMENTS:
% x   : Hurst coefficient 
% nlag: number of lags to caculate the diffrence 
% auto: the known ACF, for which we seek the Hurst coefficient    
%  
% RESULT:
% rms : the rms of the weighted differences between ideal and given ACF
%
% BACKGROUND:
%
% For a  noise series, the relationship between ACF(j) and the Hurst
% coefficient, H, is given by ACF(j) = (1/2)( |j+1|^2H -2|j|^2H + |j-1|^2H ).
%
% The Hurst estimate is reliable only for noise  series.
%===============================================================================

twoxh=min(max(2*x(1),0),1.9999);
tau=1:nfit;

%===============================================================================
% ACF(j) = (1/2)( |j+1|^2H -2|j|^2H + |j-1|^2H )
%===============================================================================
acfest=(abs(tau+1).^twoxh+(abs(tau).^twoxh)*(-2)+abs(tau-1).^twoxh)/2;
s=sum(((auto(tau+1)-acfest(tau)).^2)./tau);% weighted sum of squared differences
rms=sqrt(abs(s)/max(1,nfit));%               rms of weighted errors 

%===============================================================================
%penalty function to constrain the Hurst coefficient between 0.0001 and 0.9999
%===============================================================================
if x(1)<.0001
   rms=rms+10*abs(.0001-x(1))^.5;
end
if x(1)>.9999
   rms=rms+10*abs(.9999-x(1))^.5;
end

