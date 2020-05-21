function result=disper(ts,p);
% function result=disper(ts,p);
%
% DESCRIPTION:
% disper.m is a subfunction of FracTool, which performs  
% dispersional analysis on time series.
% ts: time series;
% p: nth power of 2, where 2^n<= length of time series;
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
% Bassingthwaighte, J.B., G.M. Raymond (1995) Evaluation of the dispersional 
% analysis method for fractal time series. Annals of Biomedical Engineering 
% 23:491-505

used=p-3;
% used: number of data points used in curve fit;
for i=1:used,
	lgsd(i)=log(stdn(ts));
	lgn(i)=log(2^i);
	tsid=(ts(1:2:length(ts)-1)+ts(2:2:length(ts)))/2;
	ts=[];
	ts=tsid;
	tsid=[];
end;
% lgsd: Natural log of SD
% lgn: Natural log of bin size
curvefit=Linreg(lgn,lgsd);
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(1)=1+curvefit(1);
result(2)=curvefit(2);
