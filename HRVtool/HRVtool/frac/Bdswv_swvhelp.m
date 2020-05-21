function result=bdswv(ts,p);
% t0=clock;
% Bridge detrended scaled window variance analysis from time series;
% ts: time series;
% p: nth power of 2, where 2^n<= length of time series;
% result(1): Hurst coefficient;
% resutt(2): correlation coefficient (r) from linear regression;
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
	sd=swvhelp(ts,n,i);   
%	counter=1;
%	sd=[];
%	for j=1:2^i:n-i,
%		tsmod=ts(j:j+2^i-1);
%		tsid=bridge(tsmod);
%		tsmod=[];
%		sd(counter)=std(tsid);
%		tsid=[];
%		counter=counter+1;
%	end;	
	lgsd(i-used1+1)=log(mean(sd));
	lgn(i-used1+1)=log(2^i);
end;
% lgsd: average standard deviation
% lgn: Natural log of bin size
curvefit=linreg(lgn,lgsd);
% curvefit: (1): slope from linear regression
%	    (2): correlation coefficient from linear regression
result(1)=curvefit(1);
result(2)=curvefit(2);
% elapsed_BDSWV_time=etime(clock,t0)