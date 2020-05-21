
% M. Noseworthy, PhD, & G.M. Wardlaw, MSc.
% A.M. Weber, MSc & M.A. Warsi, MD
% Wavelet analysis by A.M. Weber, from *insert references
% IMPORTANT- AFNI Matlab toolbox is required for THIS code.
%            get this from: http://afni.nimh.nih.gov/afni/matlab

% ========================================================================

tic
beep off;
close all;
warning off;
file = textread('List.txt', '%s', 'delimiter', '\n','whitespace', '');

for dir=1:18 %#number of subjects!
    FDdir=char(file(dir,1));
    Za=1;
    cd(FDdir)
    Bbold=strcat(FDdir, '.BRIK');
    disp(['Working on ' FDdir ' now...'])
    
    flnm=strcat(num2str(FDdir),'_FDmaps_all.mat');
    if exist(flnm, 'file')==2 %AMW
        display(' ')
        display('Loading a previously compiled full .mat file')
        eval(['load ' flnm ;]);
    else
        
        % Load BOLD data BRIK using AFNI plugin 'brikload' This will create a
        % 4D variable, BOLD, with xres x yres x #slices x time NOTE- will need
        % to keep array 'Info' as this is the header info and can be used for
        % writing BRIK files later
        
        brainbold = sprintf('%s', Bbold);   % % convert to string
        [err,BOLD,Info,ErrMessage]=BrikLoad(brainbold);
        
        % ########## TRANSPOSE DATA FOR CORRECT DISPLAY ORIENTATION
        % ############################ This rotates images 90 so they align
        % up-down so BOLD data in correct orientation: (LR,AP,SI,time)
        %newBOLD=zeros(XX,YY,ZZ,2400);
        
        BOLD=permute(BOLD,[2,1,3,4]);   % % Rearrange dimensions of array
        XX=size(BOLD,1);
        YY=size(BOLD,2);
        Zb=size(BOLD,3);
        
        BOLD=BOLD(:,:,Za:Zb,:);
        ZZ=size(BOLD,3);
        TT=size(BOLD,4);
        
        if TT==2^nextpow2(TT) %AMW_Feb2011
            powertwo=nextpow2(TT);
        else
            powertwo=nextpow2(TT)-1; % MoAl_Jan2011
        end
        bins=powertwo-2;
        
        % Calculate mean and standard deviation maps for all slices
        
        bold_mean=zeros(XX,YY,ZZ);
        bold_std=zeros(XX,YY,ZZ);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------------ Perform FD calculations on all data
        
        FDmapRDall=zeros (XX,YY,ZZ);  % % MAW_Feb2010
        FDmapPS=zeros(XX,YY,ZZ);
        FDmapWav=zeros(XX,YY,ZZ);
        RobustFDStats=zeros(XX,YY,ZZ,4); % % MAW_Feb2010 - 3rd dimension for Stats of RDfast, RDslow, RDall and FD_power (PS)
        disp('Beginning Fractal Analysis ');
        toc
        tic
        
        for SL=1:ZZ,
            BOLDtemp=squeeze(BOLD(:,:,SL,:));
            disp(['Slice # ' num2str(SL)]);
            for r=1:XX,
                for c=1:YY,
                    
                    % #### RD METHOD --------------------
                    
                    fractaldataRD=zeros(bins,2);      % % MAW_Jan2011 need to change this to (bins,2)
                    nosigcount=0;
                    % note- assumes BOLD is done with discarded acquisitions
                    % (re: steadystate).
                    databins=2^(bins+2); %Alex
                    lowb=((TT/2)-databins/2+1); %Alex
                    highb=((TT/2)+databins/2); % Alex
                    Ztemp=BOLDtemp(r,c,lowb:highb);     % % AMW_Jan2011 takes the middle pts for slice selection
                    
                    Z=(squeeze(Ztemp))'; %extract and transpose vector for pixel location 1,1.
                    
                    sumcheck=sum(Z,2);
                    
                    if sumcheck == 0,
                        % *** NO MR SIGNAL, ASSUME FD is NOISE ==1.5 *****
                        FDmapRDall (r,c,SL) = 1.5;  % % MAW_Feb2010
                        RobustFDStats(r,c,SL,:) = 0;   % % MAW_Feb2010
                        clear sumcheck;
                        nosigcount=nosigcount+1;
                    else
                        for d=0:bins, % i.e. 2^6=64       % % All data points are merged into sucessively smaller number of data points by averaging adjacent data points (i.e. they are binned)
                            m=2^d;
                            fractaldataRD(d+1,1)=m;
                            binned=reshape(Z,m,(size(Z,2)/m));      % % takes data matrix (Z) and reshapes into a new matrix depending on how many bins we need. New matrix has 'm' rows which is the number of data points to be averaged in next step
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
                    end; %end if, elseif, else logical statement
                    
                    FD_RDall = fractaldataRD(1:bins,:); % % MAW_Feb2010 - all RD time points
                    
                    % Perform fast and slow (and all) component assessment
                    
                    xmall=log10(FD_RDall(:,1)); % % MAW_Feb2010
                    ymall=log10(FD_RDall(:,2)); % % MAW_Feb2010
                    
                    % Check for NaN values using function 'isnan'
                    if sum(sum(isnan(xmall))==0 && sum(isnan(ymall))==0),
                        
                        [fitsall, stats_all]=robustfit(xmall,ymall);         % % MAW_Feb2010 - robustfit is 2X1 (y-intercept, slope)
                        FDmapRDall(r,c,SL)=1-(fitsall(2,1));    % % MAW_Feb2010 - FD(RD) = 1-slope
                        RobustFDStats(r,c,SL,3) = stats_all.p(2,1);             % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 3 for RDall)
                        
                    else
                        
                        FDmapRDall(r,c,SL)=1.5;     % % MAW_Feb2010
                        RobustFDStats(r,c,SL,:) = 0;   % % MAW_Feb2010
                    end;
                    
                    % NOW finished making: FDmapRDfast, FDmapRDslow for
                    % slice=SL % and FDmapRDall
                    
                    % ####  POWER SPEC METHOD
                    % -----------------------------------
                    fractaldataPS=zeros(bins+1,2);
                    sumcheck=sum(Z,2);
                    
                    if sumcheck == 0,
                        % *** NO MR SIGNAL, ASSUME==FD is NOISE ==1.5 *****
                        FDmapPS(r,c,SL) = 1.5;
                        RobustFDStats(r,c,SL,:) = 0;   % % MAW_Feb2010
                        clear sumcheck;
                        nosigcount=nosigcount+1;
                    else
                        Z=detrend(Z);           % % "removes the linear trend from a matrix, usually for FFT processing"
                        Yps=fft(Z, databins);       % % FFT of temporal time points
                        Pps=abs(Yps).^2;        % % Power = signal squared
                        freq=(1/(Info.TAXIS_FLOATS(2)))*((1:(databins/2))/databins);  %up to nyquist freq. %% MoAl_Jan2011 Info.TAXIS_FLOATS(2) is TR (in seconds)
                        LA=log10(Pps);
                        LA=LA(1,1:(databins/2));
                        LF=log10(freq);
                        
                        [fits, statsPS]=robustfit(LF,LA);
                        
                        %  Beta sorting for different freq spectra.  I.e one
                        %  decays as it should, one rises ... etc.
                        
                        beta=fits(2,1);
                        beta=-beta;
                        hurst=(beta+1)/2; %assuming fBn not fBm  % % actually, this fits for fGn (which is fBm)
                        D=2-hurst;	  %by selecting D=2- not 1-, scales similarly with RD.
                        %usually calculated as 2- from article.
                        FDmapPS(r,c,SL)=D;
                        RobustFDStats(r,c,SL,4) = statsPS.p(2,1);   % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 4 for PS)
                    end;
                    
                    % ####  WAVELET MLE METHOD
                    % -----------------------------------
                    
%                     cw1 = cwt(Z,1:9,'gaus2');
%                     abcw1=abs(cw1);
%                     
%                     wmat=zeros(9,512);
%                     for mgo=1:9
%                         [peaks,locs]=findpeaks(abcw1(mgo,:));
%                         sz=size(peaks);
%                         sz=sz(1,2);
%                         for go=1:sz
%                             wmat(mgo,locs(1,go))=peaks(1,go);
%                         end
%                     end
%                     
%                     ideb = 1 ; step = 1; ifin = 9;
%                     max_up=find(wmat(ideb,:));
%                     
%                     maxmat=zeros(9,512);
%                     maxmat(1,:)=wmat(1,:);
%                     
%                     for jj = ideb+step:step:ifin
%                         max_curr = find(wmat(jj,:));
%                         
%                         for k = 1:length(max_up)
%                             for l=1:length(max_curr)
%                                 if abs(max_curr(l)-max_up(k)) < 2
%                                     maxmat(jj,max_curr(l))=wmat(jj,max_curr(l));
%                                 end
%                             end
%                         end
%                         max_up = find(maxmat(jj,:));
%                     end
%                     
%                     wmat=zeros(9,512);
%                     loc=find(maxmat(1,:));
%                     for pp=1:length(loc)
%                         indeces=[1 loc(pp) maxmat(1,loc(pp))];
%                         for jj=2:9
%                             if size(indeces,1) >= jj-1
%                                 up=indeces(jj-1,2);
%                                 curr=find(maxmat(jj,:)); %
%                                 for l=1:length(curr)
%                                     if abs(curr(l)-up) == 0
%                                         newind=[jj curr(l) maxmat(jj,curr(l))];
%                                         indeces=[indeces; newind];
%                                     elseif abs(curr(l)-up) == 1
%                                         newind=[jj curr(l) maxmat(jj,curr(l))];
%                                         indeces=[indeces; newind];
%                                     end
%                                 end
%                             end
%                         end
%                         newmax=max(indeces(:,3));
%                         rows=size(indeces,1);
%                         for ii=1:rows
%                             wmat(indeces(ii,1),indeces(ii,2))=newmax;
%                         end
%                         clear indeces
%                         clear newmax
%                     end
%                     
%                     z=zeros(9,2);
%                     for pp=1:9
%                         z(pp,1)=log(sum(wmat(pp,:)));
%                         z(pp,2)=log(pp);
%                     end
%                     
%                     brob=robustfit(z(:,2),z(:,1));
%                     hurst=1+brob(2);
%                     FDmapWav(r,c,SL)=hurst;
                    
hest=wfbmesti(Z);
    hest=1+hest(3);
    FDmapWav(r,c,SL)=hest;
                    
%                     % ------------ Create plot to display signal profile [START]-----------
%                     %Create plot to display signal, RD and PS% AMW
%                     if r==int8(XX/8) && c==int8(YY/8) && SL==int8(ZZ/2)
%                         noisesig=squeeze(BOLD(r,c,SL,lowb:highb));
%                         noiseRDx=xmall;noiseRDy=ymall;
%                         noisePSf=LF;noisePSa=LA;
%                         noiseWavex=z(:,2); noiseWavey=z(:,1);
%                         noiseWaveb1=brob(1);noiseWaveb2=brob(2);
%                     end
%                     
%                     if r==int8(2*XX/5) && c==int8(2*YY/5) && SL==int8(ZZ/2)
%                         figure(1)
%                         subplot(4,2,1); plot(squeeze(BOLD(r,c,SL,lowb:highb))); title('brain signal') %*
%                         subplot(4,2,2); plot(noisesig); title('noise signal')
%                         subplot(4,2,3); plot(xmall,ymall); title('brain RD')
%                         subplot(4,2,4); plot(noiseRDx,noiseRDy); title('noise RD')
%                         subplot(4,2,5); plot(LF,LA); title('brain PS')
%                         subplot(4,2,6); plot(noisePSf,noisePSa); title('noise PS')
%                         subplot(4,2,7); scatter(z(:,2),z(:,1),'filled'); grid on; hold on
%                         plot(z(:,2),brob(1)+brob(2)*z(:,2),'g','LineWidth',2); title('brain Wave')
%                         subplot(4,2,8); scatter(noiseWavex,noiseWavey,'filled'); grid on; hold on
%                         plot(noiseWavex,noiseWaveb1+noiseWaveb2*noiseWavex,'g','LineWidth',2); title('noise Wave')
%                         eval(['print -dpsc -r1200 ' num2str(FDdir) '_signal_profile;']); %AMW
%                         close %AMW
%                     end
%                     % ------------ Create plot to display signal profile [END]-----------
                    
                    c=c+1; %column counter-----
                end;
                
                r=r+1; %row counter-----
            end;
            
            clear Ztemp;
            clear Z;
            clear sumcheck;
            
            %Graeme: Clear BOLDtemp, and re-assign for slice SL=1+SL at top of
            %loop (line 110)
            SL=SL+1;
            clear BOLDtemp;
            
        end;
        
        % ------------ Perform FD calculations [END]-----------
        
        
        % So-- save data to a matlab .mat file for later
        
        eval(['save ' num2str(FDdir) '_FDmaps_all FDdir Info XX YY ZZ FDmapRDall FDmapPS FDmapWav RobustFDStats;']);
        
    end
    
    % ----------------- Create AFNI HEAD/BRIK files [START]-----------
    % Wrote this out: AMW and MAW
    
    headers=[];
    headers=[headers 'Hurst' ];
    
    myOpt=struct();
    
    myInfo=Info;
    myInfo.DATASET_RANK(2)=1; %number of sub-bricks
    myInfo.BRICK_LABS=headers; %title for scale bar
    myInfo.TAXIS_NUMS(1) = 1; %Number of points in time (1)
    myInfo.BRICK_STATS=[0,1]; %min and max values for bricks
    myInfo.BRICK_TYPES=3; %data type: 8-byte float
    myInfo.BRICK_FLOAT_FACS=[]; % for floats, this is unscaled
    
    myOpt.Prefix=['RDall_' num2str(FDdir) ''];
    myOpt.OverWrite='y';
    Rot=zeros(64,64,ZZ);
    
    if Za==Zb       %MAW_Feb2011
        Rot(:,:,1)=2-rot90(flipud(FDmapRDall(:,:)),3);  %MAW_Feb2011
    elseif Za~=Zb
        for flp=Za:Zb
            Rot(:,:,flp)=2-rot90(flipud(FDmapRDall(:,:,flp)),3);
        end
    end
    
    [err, ErrMessage, newInfo] = WriteBrik(Rot,myInfo,myOpt);
    
    myOpt.Prefix=['PS_' num2str(FDdir) ''];
    
    if Za==Zb            %MAW_Feb2011
        Rot(:,:,1)=2-rot90(flipud(FDmapPS(:,:)),3);   %MAW_Feb2011
    elseif Za~=Zb
        for flp=Za:Zb
            Rot(:,:,flp)=2-rot90(flipud(FDmapPS(:,:,flp)),3);
        end
    end
    
    [err, ErrMessage, newInfo] = WriteBrik(Rot,myInfo,myOpt);
    
    myOpt.Prefix=['Wavf_' num2str(FDdir) ''];
    
    if Za==Zb            %MAW_Feb2011
        Rot(:,:,1)=rot90(flipud(FDmapWav(:,:)),3);   %MAW_Feb2011
    elseif Za~=Zb
        for flp=Za:Zb
            Rot(:,:,flp)=rot90(flipud(FDmapWav(:,:,flp)),3);
        end
    end
    
    [err, ErrMessage, newInfo] = WriteBrik(Rot,myInfo,myOpt);
    
    
    % ----------------- Create AFNI HEAD/BRIK files [END]-----------
    cd ../
    toc
end

