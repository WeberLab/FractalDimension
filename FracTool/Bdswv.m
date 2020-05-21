function result=bdswv(ts,p);
% function result=bdswv(ts,p);
%
% DESCRIPTION:
% bdswv.m is the subfunction of FracTool, to perform  
% bridge detrended scaled window variance analysis on time series;
% ts: fBm type time series;
% p: nth power of 2, where 2^n = length of time series;
% result(1): Hurst coefficient;
% result(2): correlation coefficient (r) from linear regression;
%
% LAST MODIFIED AT:
% 12/31 1999
%
% WRITTEN BY:
% Dr. Peter Herman (herman@elet2.sote.hu)
% Dr. Andras Eke (eke@elet2.sote.hu)
% Fractal Physiology Laboratory
% Experimental Research Department
% II. Institute of Physiology
% Semmelweis University, Budapest
% Budapest, Ulloi ut 78-A
% Hungary 1082
%
% RIGHTS AND CONDITIONS OF USE:
% This is freeware. Please place the following reference in any publication for which 
% this software or any derivates of it is used and send one reprint to Dr. Eke 
% at the address given above:
%
% A. Eke, P. Hermán, J. B. Bassingthwaighte, G. M. Raymond, D. B. Percival, 
% M. Cannon, I. Balla, and C. Ikrényi. Physiological time series: distinguishing
% fractal noises from motions. Pflügers Archiv European Journal of Physiology, 
% 439(4):403-415, 2000
%
% This software is permitted to be distributed only as a complete compressed package
% (see below). The README file must be included.
%
% ACKNOWLEDGEMENTS:
% FracTool was developed with support from OTKA Grants I/3 2040, T 016953, 
% and NIH grant TW00442. 
%
% CONTENTS OF FRACTOOL PACKAGE:
% README.txt
% Bridge.m
% Disper.m
% Bdswv.m
% Stdn.m
% Linreg.m
% Spec.m
% Window.m
% Fractool.m
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
		tsid=Bridge(tsmod);
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
curvefit=Linreg(lgn,lgsd);
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(1)=curvefit(1);
result(2)=curvefit(2);