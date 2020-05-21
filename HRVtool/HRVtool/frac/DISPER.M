function result=disper(ts,p);
% function result=disper(ts,p);
%
% DESCRIPTION:
% disper.m is the subfunction of FractX, to evaluate  
% dispersional analysis of time series.
% ts: time series;
% p: nth power of 2, where 2^n<= length of time series;
% result(1): Hurst coefficient;
% result(2): correlation coefficient (r) from linear regression;
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
% Bassingthwaighte, J.B., G.M. Raymond (1995) Evaluation of the dispersional 
% analysis method for fractal time series. Annals of Biomedical Engineering 
% 23:491-505

used=p-3;
% used: number of data points used in curve fit;
for i=1:used,
	%lgsd(i)=log(stdn(ts));
   lgsd(i)=log(std(ts));
   lgn(i)=log(2^i);
	tsid=(ts(1:2:length(ts)-1)+ts(2:2:length(ts)))/2;
	ts=[];
	ts=tsid;
	tsid=[];
end;
% lgsd: Natural log of SD
% lgn: Natural log of bin size
curvefit=linreg(lgn,lgsd);
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(1)=1+curvefit(1);
result(2)=curvefit(2);
