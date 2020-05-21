function result=spec(ts,p,f,avg);
% function result=spec(ts,p,f,avg);
% Spectral analysis of time series;
% ts: time series;
% p: nth power of 2, where 2^n<= length of time series;
% f: mintaveteli frekvencia;
% avg: if avg=0, then spectra all, if avg=1, then spectra avg;
% result(1): Hurst coefficient;
% result(2): correlation coefficient (r) from linear regression;
% result(3): beta from Fourier analysis;
% result(4): Hurst coefficient (not used high frequencies);
% result(5): correlation coefficient (r) from linear regression (not used high frequencies);
% result(6): beta from Fourier analysis (not used high frequencies);

%t0=clock;
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
curvefit=linreg(lgfrek,lgPyy);
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
   curvefit=linreg(lgfrek(1:length(lgfrek)-2),lgPyy(1:length(lgPyy)-2));
else
   curvefit=linreg(lgfrek(1:floor(length(lgfrek)/4)),lgPyy(1:floor(length(lgPyy)/4)));
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