function fbn=dhm(nlog,sed,h);
% function fbn=dhm(nlog,sed,h);
% Fractional Gaussian process method using Davies-Harte method.
% fbn: fractional Brownian noise;
% nlog: the length of time series=2^nlog;
% sed: seed for the random number generator;
% h: Hurst coefficient;
% sd: initial standard deviation;
t0=clock;
n=2^nlog;
sd=1;

%Generate autocovariance sequence;
var=sd*sd;
M=2*n;
r=1:n;
hurst=2*h;
s(2:n+1)=var/2.*(abs(r+1).^hurst-2*abs(r).^hurst+abs(r-1).^hurst);
s(1)=var;
clear r;

%Compute SV(j) from the discrete Fourier transform of s(j);
s(n+2:M)=s(n:-1:2);
% s az SV helyebe lep memoria takarekossagbol
s=fft(s);

s=real(s);
%for i=1:M,
%	if s(i)<0, s(i)=0; end;
%end;
s(find(s<0))=0;
s=s.^0.5;

%Generate the random numbers with zero mean and unit variance;
randn('seed',sed);
w=randn(1,M);


%Calculate V(j);
sqroh=0.5^0.5;
V(1)=s(1)*w(1);
V(2:n)=sqroh.*s(2:n).*(w(2:2:M-2)+w(3:2:M-1).*j);
V(n+1)=s(n+1)*w(M);
clear s w;
V(n+2:M)=conj(V(n-1:-1:1));

%Compute fbn;

V=fft(V);
sqM=M^0.5;

fbnr=real(V(1:n))'./sqM;
fbni=imag(V(1:n))'./sqM;
clear V;
fbn=fbnr+fbni;

elapsed_FGP_time=etime(clock,t0)
