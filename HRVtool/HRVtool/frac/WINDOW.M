function result=window(ts);
% function result=window(ts);
% Windowing for spectral analysis of time series;
% ts: time series;
%t0=clock;
n=length(ts);
w=1-(2*(1:n)/(n+1)-1).^2;
w=w';
result=ts.*w;
%elapsed_windowing_time=etime(clock,t0)