function ps=cgsa(ts,H);
%function cgsa(ts,H);
%ts: timeseries
%H:apriori knowledge of Hurst-coefficient
%ps: power spectra of 1/f noise-free time series

l=length(ts);
xx=zeros(2*l,1);
xx(1:2:end-1)=ts;
xx(2:2:end)=ts;
tsf=xx(1:floor(end/2));
tsf=tsf./2^H;
fx=fft(ts);
fxf=fft(tsf);
ps=fx.*conj(fx)-abs(fx.*conj(fxf));