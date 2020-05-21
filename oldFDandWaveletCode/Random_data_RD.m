a=-pi:0.0399:25*pi;         
A=sin(a(1:2048));           % sin wave

b=-2000*pi:1:2000*pi;
B=sin(b(1:2048));           % different sin wave

C = rand(1,2048);           % random numbers

D = ones(1,2048);           % flat line

E = (1:2048);               % line with gradient = 1

F = [1:0.5:255,256:2:1620,1620.5:0.5:2048];  % non-straight line with gradient approx = 1

G = [1:512,(200*ones(1,1023)),512:1024];   % very messy line with gradient approx = 0.5

Z = A;          % Change Z = A or B or C etc depending on what data you want to run

plot (Z)


for d=0:9, % i.e. 2^9=512
				m=2^d;
				fractaldataRD(d+1,1)=m;
				binned=reshape(Z,m,(size(Z,2)/m));
				sigbin=mean(binned,1); %binning consecutive values

				SD_bin=std(sigbin,1); %normalise by N not N-1
				meanbin=mean(sigbin);
						
				RD=SD_bin/meanbin; % % for current bin size
			 	fractaldataRD(d+1,2) = RD;

				clear binned;
				clear sigbin;
				clear meanbin;
				clear SD_bin;
				clear RD;
			
			    d=d+1;	
end


FD_RDfast = fractaldataRD(1:4,:); % 1st 4 points used in FDfast fit
			FD_RDslow = fractaldataRD(5:9,:); % 1ast 5 points used in FDslow fit
            FD_RDall = fractaldataRD; % % MAW_Feb2010 - all RD time points
            

            % Perform fast and slow (and all) component assessment
            
            xmfast=log10(FD_RDfast(:,1));
			ymfast=log10(FD_RDfast(:,2));
			xmslow=log10(FD_RDslow(:,1));
			ymslow=log10(FD_RDslow(:,2));
            xmall=log10(FD_RDall(:,1)); % % MAW_Feb2010
            ymall=log10(FD_RDall(:,2)); % % MAW_Feb2010
            
fit=robustfit(xmall,ymall);
alpha=fit(2,1);
RD=1-alpha;
H_RD=alpha-1;

%% Power spetrum

TR=0.25;

Yps=fft(Z, 2048);       % % FFT of temporal time points
			Pps=abs(Yps).^2;        % % Power = signal squared
			freq=(1/(TR))*((1:1024)/2048);  %up to nyquist freq.
			LA=log10(Pps);
			LA=LA(1,1:1024);
			LF=log10(freq);

			[fits, statsPS]=robustfit(LF,LA);

%  Beta sorting for different freq spectra.  I.e one decays as it should, one rises ... etc.

			beta=fits(2,1);
			beta=-beta;
			hurst=(beta+1)/2; %assuming fBn not fBm  % % actually, this fits for fGn (which is fBm)
			PS=2-hurst;	  %by selecting D=2- not 1-, scales similarly with RD.
				          %usually calculated as 2- from article.
            H_PS=hurst;
                          
%%

PS
RD
