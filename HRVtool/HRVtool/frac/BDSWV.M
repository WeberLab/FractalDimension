function result=bdswv(ts,p);
% function result=bdswv(ts,p);
%
% DESCRIPTION:
% bdswv.m is the subfunction of FractX, to evaluate  
% bridge detrended scaled window variance analysis from time series;
% ts: fBm type time series;
% p: nth power of 2, where 2^n = length of time series;
% result(1): Hurst coefficient;
% result(2): correlation coefficient (r) from linear regression;
%
% FractX is a mosaic name composed from Single (Fract)al Inde(X), 
% a new fractal measure.
%
% LAST MODIFIED AT:
% 02/09 1998
%
% WRITTEN BY:
% Dr. Peter Herman
% Dr. Andras Eke
% Fractal Physiology Lab
% Experimental Research Department
% II. Institute of Physiology
% Budapest, Ulloi ut 78-A
% Hungary 1082
%
% RIGHTS AND ACKNOWLEDGEMENTS:
% This software was developed with support from OTKA Grants I/3 2040, T 016953, 
% and NIH grant TW00442. Please cite these grants in any publication for which 
% this software or any derivates of it is used and send one reprint to Dr. Eke
% at the address given above.
%
% REFERENCE:
% Cannon, M.J., D.B. Percival, D.C. Caccia, G.M. Raymond, J.B. Bassingthwaighte
% (1997) Evaluating scaled window variance methods for estimating the Hurst
% coefficient of time series. Physica A 241:606-626

ignored=[6 1 0; 7 1 0; 8 1 0; 9 2 2; 10 2 3; 11 2 4; 12 2 4; 13 2 5; 14 3 6; 15 3 7; 16 3 7; 17 3 7; 18 3 7];
% ignored: length, smallest bin size, largest bin size;
i=1;
while p>=ignored(i,1),
	smallest=ignored(i,2);
	largest=ignored(i,3);
	i=i+1;
end;
% smallest: ignored these smallest bin size;
% largest: ignored these largest bin size;
used1=smallest+1;
used2=p-largest;
n=length(ts);
lgsd=zeros(1,used2-used1+1);
lgn=zeros(1,used2-used1+1);
for i=used1:used2,
	counter=1;
	sd=[];
	for j=1:2^i:n-i,
		tsmod=ts(j:j+2^i-1);
		tsid=bridge(tsmod);
		tsmod=[];
		sd(counter)=std(tsid);
		tsid=[];
		counter=counter+1;
	end;	
	lgsd(i-used1+1)=log(mean(sd));
	lgn(i-used1+1)=log(2^i);
end;
% lgsd: Natural log of average standard deviation
% lgn: Natural log of bin size
curvefit=linreg(lgn,lgsd);
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(1)=curvefit(1);
result(2)=curvefit(2);