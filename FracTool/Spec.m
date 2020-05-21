function result=spec(ts,p,avg);
% function result=spec(ts,p,f,avg);
% Spectral analysis of time series;
% ts: time series;
% p: nth power of 2, where 2^n<= length of time series;
% avg: if avg=0, then spectra all, if avg=1, then spectra avg;
% result(1): Hurst coefficient;
% result(2): correlation coefficient (r) from linear regression;
% result(3): beta from Fourier analysis;
% result(4): Hurst coefficient (high frequencies excluded);
% result(5): correlation coefficient (r) from linear regression (high frequencies excluded);
% result(6): beta from Fourier analysis (high frequencies excluded);
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

%t0=clock;
f=1;
n=2^p;
Y = fft(ts);
Pyy = Y.*conj(Y)/n;
Pyy=Pyy(1:n/2);
Y=[];y=[];
frek = f/n*(1:n/2);
frek=frek';
if avg==1,
	for i=1:(p-1),
		ff(i)=mean(frek(2^(i-1):2^i-1));
		pp(i)=mean(Pyy(2^(i-1):2^i-1));
	end;
	Pyy=[];
	Pyy=pp(2:length(pp));
	frek=[];
	frek=ff(2:length(pp));
end;
k=find(Pyy==0);
if isempty(k)==0, 
   for kk=1:length(k)
      Pyy(k(kk))=0.00001;
   end;
end;
lgPyy=log(Pyy);
lgfrek=log(frek);
% lgfrek: Natural log of frequency
% lgPyy: Natural log of Power spectrum
curvefit=Linreg(lgfrek,lgPyy);
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(3)=curvefit(1)*(-1);
result(2)=curvefit(2);
if result(3)<1,
	result(1)=(result(3)+1)/2;
end;
if result(3)>1,
	result(1)=(result(3)-1)/2;
end;
if result(3)==1,
	result(1)=inf;
end;
if avg==1,
   curvefit=Linreg(lgfrek(1:length(lgfrek)-2),lgPyy(1:length(lgPyy)-2));
else
   curvefit=Linreg(lgfrek(1:floor(length(lgfrek)/4)),lgPyy(1:floor(length(lgPyy)/4)));
end;
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(6)=curvefit(1)*(-1);
result(5)=curvefit(2);
if result(6)<1,
	result(4)=(result(6)+1)/2;
end;
if result(6)>1,
	result(4)=(result(6)-1)/2;
end;
if result(6)==1,
	result(4)=inf;
end;
%elapsed_spect_time=etime(clock,t0)