
hestlst=zeros(100,3);
rdlst=zeros(100,3);
pslst=zeros(100,3);
ekelst=zeros(100,3);

jake=1
%for jake=1:100
    
    rndn=rand(1);
    %FBW=wfbm(rndn,512);
    %FBW = ffgn(1,rndn,1,512,0);
    otf=1:1:512;
    FBW=sin(otf);
    noise=awgn(FBW,1);
    FBW=noise;
    
    %RD
    fractaldataRD=zeros(9,2);
    for d=0:9,
        m=2^d;
        fractaldataRD(d+1,1)=m;
        binned=reshape(FBW,m,(size(FBW,2)/m));      % % takes data matrix (Z) and reshapes into a new matrix depending on how many bins we need. New matrix has 'm' rows which is the number of data points to be averaged in next step
        sigbin=mean(binned,1); %binning consecutive values      % % the reshaped columns are averaged to give us our new bins (thus a two step process of going from data (Z) to a series of averaged adjacent data points)
        
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
    end; % endfractaldim calc loop for pixel timecourse
    FD_RDall = fractaldataRD(1:9,:);
    xmall=log10(FD_RDall(:,1)); % % MAW_Feb2010
    ymall=log10(FD_RDall(:,2)); % % MAW_Feb2010
    [fitsall, stats_all]=robustfit(xmall,ymall);         % % MAW_Feb2010 - robustfit is 2X1 (y-intercept, slope)
    FDmapRDall=2-(1-(fitsall(2,1)));
    
    rdlst(jake,:)=[FDmapRDall];
    
    
    %PS
    fractaldataPS=zeros(9+1,2);
    Z=detrend(FBW);           % % "removes the linear trend from a matrix, usually for FFT processing"
    Yps=fft(Z, 512);       % % FFT of temporal time points
    Pps=abs(Yps).^2;        % % Power = signal squared
    freq=(1/2)*((1:(512/2))/512);  %up to nyquist freq. %% MoAl_Jan2011 Info.TAXIS_FLOATS(2) is TR (in seconds)
    LA=log10(Pps);
    LA=LA(1,1:(512/2));
    LF=log10(freq);
    
    [fits, statsPS]=robustfit(LF,LA);
    
    beta=fits(2,1);
    beta=-beta;
    hurstPS=((beta+1)/2); %assuming fBn not fBm  % % actually, this fits for fGn (which is fBm)

    pslst(jake,:)=[hurstPS];
    
    %wfbmesti
    
    hest=wfbmesti(FBW);
    hest=1+hest(3);
    
    hestlst(jake,:)=[hest];
    
    %eke
    Z=FBW;
    orig_signal=FBW;
    time_points=512;
    bins=7;
    meanZ=mean(FBW);
    j=1:time_points;
    
    % End match (subtract line connecting first and last points from data)
    y1 = Z(1); y2 = Z(time_points);
    slope = (y2-y1)/(time_points-1);
    y_int = y2 - slope*(time_points);      % first time point is 1, or 0?
    %for k=1:n
    E = slope*j + y_int;
    Z = Z - E;
    
    % Multiply each new value by parabolic window
    W = 1 - (2*j/(time_points+1)-1).^2;    % parabolic window
    Z = Z.*W;
    
    % Power spectrum
    fs = 1.0/2;
    range = ceil((time_points+1)/2); % half, including f_nyquist
    
    fftSignal = fft(Z,time_points);
    fftSignal = fftSignal(1:range);  % 1st half of fft since it's symmetric
    
    PSD = abs(fftSignal).^2/time_points;
    %               if rem(n,2)
    %                  PSD(2:end) = PSD(2:end)*2;    % double PSD amplitude at non-unique points
    %               else
    %                  PSD(2:end-1) = PSD(2:end-1)*2;
    %               end
    
    freq = fs*(0:range-1)/time_points;         % vector for frequency axis
    
    % Evaluate only low frequencies (1/8 < fs < 1/2) if requested
    
    PSDpart=PSD(2:time_points/8);
    freqPart=freq(2:time_points/8);
    
    % switching back to gw code
    
    LA = log10(PSDpart);
    LF = log10(freqPart);
    
    [fits, statsPS]=robustfit(LF,LA);
    beta=-fits(2,1);
    
    % now we have our low freq end matched beta value
    
    %temp to see if this thing works
    hurst=(beta+1)/2; %assuming fBn not fBm  % % actually, this fits for fGn (which is fBm)
    D=2-hurst;	  %by selecting D=2- not 1-, scales similarly with RD.
    %usually calculated as 2- from article.
    
    ssc_result=0;
    
    if beta >= 0.38 && beta <=1.04
        %do SSC
        
        Z=orig_signal;
        % cumulative sum of signal
        for j = 2:time_points
            Z(j) = Z(j) + Z(j-1);
        end
        
        % bridge detrending
        j = 1:time_points;
        y1 = Z(1); y2 = Z(time_points);
        slope = (y2-y1)/(time_points-1);
        %y_int = signal(1);
        E = slope*j + y1;
        Z = Z - E;
        
        % SWV
        SWV = zeros(bins,1);
        mVec = zeros(bins,1);
        for d = 0:bins, % i.e. 2^9=512
            m = 2^d;
            mVec(d+1) = m;
            binned = reshape(Z, [m, (length(Z)/m)]);
            stdBin = std(binned(1:end,:));
            SWV(d+1) = mean(stdBin);
        end
        
        SWVlog = log10(SWV); % would keep pts 3,4,5,6,7 but robust fit will take care of that
        mlog = log10(mVec);
        
        [fits, statsPS]=robustfit(mlog,SWVlog);
        ssc_h=fits(2,1);
        hurst=ssc_h;
        
        if ssc_h >0.9 % actually hurst 0.9 +-0.1 is unclasifiable
            ssc_result=1; % i.e signal should be analysed as fBm
        else
            ssc_result=0; % i.e signal should be analysed as fGn
        end
    end
    
    
    if beta<0.38 || ssc_result==0
        hurst=(beta+1)/2;
        D=2-hurst;

        % RD is as calculated above
    elseif beta>1.04 || ssc_result==1
        hurst=(beta-1)/2;
        D=2-hurst;
        
        %do SWV to get new RD value
        % like SSC without cumalative sum
        
        Z=orig_signal;
        % bridge detrending
        j = 1:time_points;
        y1 = Z(1); y2 = Z(time_points);
        slope = (y2-y1)/(time_points-1);
        %y_int = signal(1);
        E = slope*j + y1;
        Z = Z - E;
        
        % SWV
        SWV = zeros(bins,1);
        mVec = zeros(bins,1);
        for d = 0:bins, % i.e. 2^9=512
            m = 2^d;
            mVec(d+1) = m;
            binned = reshape(Z, [m, (length(Z)/m)]);
            stdBin = std(binned(1:end,:));
            SWV(d+1) = mean(stdBin);
        end
        
        SWVlog = log10(SWV); % would keep pts 3,4,5,6,7 but robust fit will take care of that
        mlog = log10(mVec);
        
        [fits, statsPS]=robustfit(mlog,SWVlog);
        H_rd=fits(2,1);       
    end

    ekelst(jake,:)=[hurst];
    
    clearvars -except rdlst pslst hestlst ekelst
%end 

%stats=[mean(rdlst(:,3)) std(rdlst(:,3)); mean(pslst(:,3)) std(pslst(:,3)); mean(hestlst(:,3)) std(hestlst(:,3)); mean(ekelst(:,3)) std(ekelst(:,3))];
    
