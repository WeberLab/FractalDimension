% for me to do: automation and ROI with masking for alex to do: fix
% component section and ROI other than masking


% M. Noseworthy, PhD, & G.M. Wardlaw, MSc. IMPORTANT- AFNI Matlab toolbox
% is required for THIS code.
%            get this from: http://afni.nimh.nih.gov/afni/matlab This also
%            includes info on how HEAD/BRIK load/save works

%  HERE's what this program assumes you have already done: 1.  The time
%  series data has been renamed to I.#.# scheme using AFNIrename 2.  You
%  have used AFNI "to3d" to save BOLD data as HEAD and BRIK files
% e.g. to3d -epan -time:tz 2400 3 250 altplus -save_outliers bad.txt
% Ser3/I.0*.*

% NOTE this uses MATLAB AFNI toolbox functions 'BrikLoad' and 'WriteBrik';

% ========================================================================
% MDN Sept 9/09  Used permute command to rotate matrix dimensions (thx Arv)
% MDN Sept 9/09  Added functionality to work with any sie BOLD matrices

% ========================================================================

% % MAW 14 Dec 2009 V4 changes MRSy to + from - since the input Y %
% coordinates now have the correct sign.

% % MAW_Feb2010 % now takes data points from the MIDDLE of the data set %
% Also, added stats for RD and PS fits, added a new matrix of these stats %
% Added a new RD fit, with all data (not assume two trends) % changed the
% way the figures are produced % changed the stats to include the number of
% significant fits

% % MoAl_Jan2011 % Major overhall of program by AMW and MAW to allow for
% various input % data types and various analyses. Allows for whole brain
% analysis % (TR=2000 ms) as well as traditional 3 slice (TR=250 ms).

% % Note: comments added by MAW usually begin with '% %'


% % MAW_Feb2011
function FDall_MoAl_MAWV05_AMWV05(auto) % note, this line must ALWAYS match the filename
% Runs this program as a function asking for automation
% run this by calling FDall_MoAl_MAWV04_AMWV05(y/n)

if exist('auto','var')==0     % % MAW_Feb2011
    automation='n';
else
    automation=auto;
end

if automation=='y'      % % MAW_Feb2011
    global FDdir
    global Za
    global Zb
    global doSlice
    global doComponents
    global doFD
    global Bbold
    global bins
    global gofig
    global ROIyes
    global ROItype
    global roi_number
    global ROIname
    global maskfile
    global rn
    global ln
    global an
    global pn
    global in
    global sn
    global SLnum
end

ROI_size=0;

tic
beep off;
close all;
warning off;

% automation='n'; %% MAW_Jan2011 change this to 'y' if running in a script

if automation=='n'  %% MoAl_Jan2011
    disp('This program is an automated multislice FD analysis'); % saving as Matlab HEAD/BRIK files');
    disp(' ');
    disp('It has 3 parts:');
    disp('   1) Creates .mat datafiles of various FD parameters');
    disp('   2) Saves high res postscript files of FD parametric maps');
    disp('   3) Performs optional ROI analysis');
    disp('It can work on previously saved .mat files of FD data as well');
    disp(' ');
    disp('TYPE ANY KEY to continue to select directory with your data (HEAD and BRIK files)');   %% MoAl_Jan2011
    disp(' ')
    pause
    cd(uigetdir)
    ls
    FDdir=input('Enter output .mat file prefix (may use exam ID): ','s');
    
    corin=0;
    while (corin==0); %AMW, idiot proof
        disp(' ');
        disp('Would you like to analyse All slices, a Range of slices')
        disp('or a Single slice (a/r/s)?')
        doSlice=input('OR type afni to enter afni and then decide: ','s');   %% MoAl_Jan2011
        if doSlice=='s'                                            %% MoAl_Jan2011
            Za=input('Which slice number? ');
            Zb=Za;
            corin=1;
        elseif doSlice=='r'
            Za=input('First slice? ');
            Zb=input('Last slice? ');
            corin=1;
        elseif doSlice=='a'
            Za=1;
            corin=1;
        elseif doSlice=='afni'
            !afni &
        else
            disp('You did not enter a, r or s (or afni)');
        end
    end
    
    corin=0;
    while (corin==0); %AMW, idiot proof
        disp(' ');
        doComponents=input('Would you like to analyse Fast and Slow components (y/n)? ','s');
        if doComponents=='y' || doComponents=='n'
            corin=1;
        else
            disp('You did not enter y or n');
        end
    end
    
elseif automation=='y'      %% MAW_Feb2011
    disp('Running in automation mode');
    
else
    disp('You did not choose y/n for automation'); %AMW, Mo, you might want to use a WHILE statement for this
end

if automation=='n'  % MAW
    corin=0;
    while (corin==0) %AMW
        disp(' ');
        disp('There may already be some FD maps in this directory');
        doFD=input('Is there a .mat FDmap file already saved that you would like to use (y/n)? ','s');
        if doFD=='y' || doFD=='n'
            corin=1;
        else
            display('You are doing it wrong')
        end
    end
end

if doFD=='n',
    if automation=='n';
        display(' ')
        display('TYPE ANY KEY to continue to select the file you would like to analyze (BRIK file)')
        display('note: this may take a while, and MATLAB may continue to display "Pause" even though its actually Busy')
        pause
        Bbold=uigetfile('*.BRIK', 'Select the BRIK file you would like to analyze'); %AMW
    end
    
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
    
    if doSlice=='a'
        Zb=size(BOLD,3);
    end
    
    BOLD=BOLD(:,:,Za:Zb,:);
    ZZ=size(BOLD,3);
    TT=size(BOLD,4);
    
    if TT==2^nextpow2(TT) %AMW_Feb2011
        powertwo=nextpow2(TT);
    else
    powertwo=nextpow2(TT)-1; % MoAl_Jan2011
    end
    
    disp(' ')
    disp(['Your data has a max set of ' num2str(powertwo-2) ' bins.']); % MoAl_Jan2011
    
    if automation=='n'
        corin=0;
        while (corin==0) %AMW
            bins=input('How many bins would you like to analyze? '); % MoAl_Jan2011
            if bins <= powertwo-2
                corin=1;
            else
                display ('Try again bozo')
            end
        end
    elseif automation=='y'  % % MAW_Feb2011
        if bins==99 % MAW_Feb2011 if want to use max number of bins
            bins=powertwo-2;
        end
    end
    
    % Calculate mean and standard deviation maps for all 3 slices
    
    bold_mean=zeros(XX,YY,ZZ);
    bold_std=zeros(XX,YY,ZZ);
    
    for cv=1:ZZ,
        bold_mean(:,:,cv)=mean(squeeze(BOLD(:,:,cv,:)),3);     % % squeeze: Remove singleton dimensions
    end;
    clear cv;
    
    for cv=1:ZZ,
        bold_std(:,:,cv)=std(squeeze(BOLD(:,:,cv,:)),1,3);
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------------ Perform FD calculations on all data
    % [START]------------------ Here we are fitting 2 components- the fast
    % and slow FDs.
    
    %Assigning zero filled matrices for FD RD fast, slow and FD PS
    FDmapRDall=zeros (XX,YY,ZZ);  % % MAW_Feb2010
    
    if doComponents=='y'            % MoAl_Jan2011
        FDmapRDfast=zeros(XX,YY,ZZ);
        FDmapRDslow=zeros(XX,YY,ZZ);
        midbins=round(bins/2);
    end
    
    FDmapPS=zeros(XX,YY,ZZ);
    RobustFDStats=zeros(XX,YY,ZZ,4); % % MAW_Feb2010 - 3rd dimension for Stats of RDfast, RDslow, RDall and FD_power (PS)
    anat_size=[512 512];  %Change this when and if anat briks are loaded in for interps/overlays later
    interpsize = anat_size(1);
    disp('Beginning Fractal Analysis ');
    toc
    tic
    
    for SL=1:ZZ,
        BOLDtemp=squeeze(BOLD(:,:,SL,:));
        
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
                    FDmapRDfast(r,c,SL) = 1.5;
                    if doComponents=='y'
                        FDmapRDslow(r,c,SL) = 1.5;
                        FDmapRDall (r,c,SL) = 1.5;  % % MAW_Feb2010
                    end
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
                
                if doComponents=='y'
                    FD_RDfast = fractaldataRD(1:(midbins-1),:); % 1st 4 points used in FDfast fit     % % more data points, smaller bins ('m' is low)
                    FD_RDslow = fractaldataRD(midbins:bins,:); % 1ast 5 points used in FDslow fit    % % fewer data points because large bins have been averaged
                end
                FD_RDall = fractaldataRD(1:bins,:); % % MAW_Feb2010 - all RD time points
                
                % Perform fast and slow (and all) component assessment
                if doComponents=='y'
                    xmfast=log10(FD_RDfast(:,1));   % % first column is 'm' which was number in each bin
                    ymfast=log10(FD_RDfast(:,2));   % % second column was the RD
                    xmslow=log10(FD_RDslow(:,1));
                    ymslow=log10(FD_RDslow(:,2));
                end
                xmall=log10(FD_RDall(:,1)); % % MAW_Feb2010
                ymall=log10(FD_RDall(:,2)); % % MAW_Feb2010
          
                % Check for NaN values using function 'isnan'
                if doComponents=='y'
                    if sum(isnan(xmfast))==0 && sum(isnan(ymfast))==0 && sum(isnan(xmslow))==0 && sum(isnan(ymslow))==0,
                        [fitsfast, stats_fast]=robustfit(xmfast,ymfast);      % % MAW_Feb2010 added stats meassure to each fit
                        FDmapRDfast(r,c,SL)=1-(fitsfast(2,1));
                        RobustFDStats(r,c,SL,1) = stats_fast.p(2,1);            % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 1 for RDfast)
                        % % Note that stats.p gives P value for the
                        % y-intercept and for the slope (2X1 matrix)
                        [fitsslow, stats_slow]=robustfit(xmslow,ymslow);
                        FDmapRDslow(r,c,SL)=1-(fitsslow(2,1));
                        RobustFDStats(r,c,SL,2) = stats_slow.p(2,1);            % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 2 for RDslow)
                        
                    else
                        FDmapRDfast(r,c,SL)=1.5;
                        FDmapRDslow(r,c,SL)=1.5;
                        
                        RobustFDStats(r,c,SL,:) = 0;   % % MAW_Feb2010
                    end;
                end
                
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
                
                % ------------ Create plot to display signal profile [START]-----------
                %Create plot to display signal, RD and PS% AMW
                if r==int8(XX/8) && c==int8(YY/8) && SL==int8(ZZ/2)
                    noisesig=squeeze(BOLD(r,c,SL,lowb:highb));
                    noiseRDx=xmall;noiseRDy=ymall;
                    noisePSf=LF;noisePSa=LA;
                end
                
                if r==int8(2*XX/5) && c==int8(2*YY/5) && SL==int8(ZZ/2)
                    figure(1)
                    subplot(3,2,1); plot(squeeze(BOLD(r,c,SL,lowb:highb))); title('brain signal') %*
                    subplot(3,2,2); plot(noisesig); title('noise signal')
                    subplot(3,2,3); plot(xmall,ymall); title('brain RD')
                    subplot(3,2,4); plot(noiseRDx,noiseRDy); title('noise RD')
                    subplot(3,2,5); plot(LF,LA); title('brain PS')
                    subplot(3,2,6); plot(noisePSf,noisePSa); title('noise PS')
                    eval(['print -dpsc -r1200 ' num2str(FDdir) '_signal_profile;']); %AMW
                    close %AMW
                end
                % ------------ Create plot to display signal profile [END]-----------
                
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
    
    disp('Starting image interpolation');
    toc
    tic
    
    % DO These interpolations After Loop!!
    if doComponents=='y'
        
        FDmapRDfast_interp=imresize(FDmapRDfast,[interpsize interpsize],'bicubic');
        FDmapRDslow_interp=imresize(FDmapRDslow,[interpsize interpsize],'bicubic');
    end
    
    FDmapRDall_interp=imresize(FDmapRDall,[interpsize interpsize],'bicubic');
    FDmapPS_interp=imresize(FDmapPS,[interpsize interpsize],'bicubic');
    
    
    % ------------ Perform FD calculations [END]-----------
    
    
    % So-- save data to a matlab .mat file for later
    if doComponents=='y' %AMW changed naming, to include what option the user used (bins, components, slices)
        if doSlice=='a' && bins==powertwo
            eval(['save ' num2str(FDdir) '_FDmaps_all_Components FDdir Info XX YY ZZ bold_mean bold_std FDmapRDslow FDmapRDfast FDmapRDall FDmapPS RobustFDStats FDmapRDslow_interp FDmapRDfast_interp FDmapRDall_interp FDmapPS_interp;']);
        else
            eval(['save ' num2str(FDdir) '_FDmaps_Slices_' num2str(Za) '-' num2str(Zb) '_Bins_' num2str(bins) '_Components FDdir Info XX YY ZZ bold_mean bold_std FDmapRDslow FDmapRDfast FDmapRDall FDmapPS RobustFDStats FDmapRDslow_interp FDmapRDfast_interp FDmapRDall_interp FDmapPS_interp;']);
        end
    else
        if doSlice=='a' && bins==powertwo
            eval(['save ' num2str(FDdir) '_FDmaps_all FDdir Info XX YY ZZ bold_mean bold_std FDmapRDall FDmapPS RobustFDStats FDmapRDall_interp FDmapPS_interp;']);
        else
            eval(['save ' num2str(FDdir) '_FDmaps_Slices_' num2str(Za) '-' num2str(Zb) '_Bins_' num2str(bins) ' FDdir Info XX YY ZZ bold_mean bold_std FDmapRDall FDmapPS RobustFDStats FDmapRDall_interp FDmapPS_interp;']);
        end
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
    
    %need to add for components...laaazy%
    
    % ----------------- Create AFNI HEAD/BRIK files [END]-----------
    
else %AMW option to load desired .mat file. What about automation? Easy for previously compiled full (all slices, all bins) .mat files...
    if doComponents=='y' %AMW
        flnm=strcat(num2str(FDdir),'_FDmaps_all_Components.mat');
    else
        flnm=strcat(num2str(FDdir),'_FDmaps_all.mat');
    end
    
    if exist(flnm, 'file')==2 %AMW
        display(' ')
        display('Loading a previously compiled full .mat file')
        eval(['load ' flnm ;]);
    else
        display(' ')
        display('TYPE ANY KEY to select the .mat file you would like to use')
        pause
        ldmat=uigetfile('*.mat', 'Select the .mat file you would like to use'); %AMW
        eval(['load ' ldmat ;]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parts below are optional and prompt the user for decision as to
% whether they are to be done or not. ----------------- OPTIONAL Save
% Figures [START]-----------

if automation=='n'
    corin=0;
    while (corin==0) %AMW
        display(' ')
        gofig=input('Do you want to make high res ps image files (y/n) ? ','s');
        if gofig=='y' || gofig=='n'
            corin=1;
        else
            display('well, do you or not?')
        end
    end
else
    gofig='y';
end

if gofig == 'y',
    if doComponents=='y'
        
        for a=1:ZZ;  % % MoAl_Jan2011 %AMW
            
            % Create parametric maps from each slice-------------
            figure(a); %AMW
            subplot(5,2,1);         % % Sub divide figure in a 1 X 2 'sub plots' and place this one in box 1
            imagesc(bold_mean(:,:,a)); axis('image'); colormap('gray'); axis off;
            title('BOLD MEAN');
            subplot(5,2,2);
            imagesc(bold_std(:,:,a)); axis('image'); colormap('gray'); axis off;
            title('BOLD STDev');
            subplot(5,2,3);         % % Note, fills sub plots across row then down to next column
            imagesc(FDmapRDfast_interp(:,:,a));colormap('jet');
            colorbar; axis('image'); axis off;
            title('FD Map: Relative Dispersion:FAST --Interpolated');
            subplot(5,2,4);
            imagesc(FDmapRDslow_interp(:,:,a));colormap('jet');
            colorbar; axis('image'); axis off;
            title('FD Map: Relative Dispersion:SLOW --Interpolated');
            subplot(5,2,5);
            imagesc(FDmapRDall_interp(:,:,a));colormap('jet');
            colorbar; axis('image'); axis off;
            title('FD Map: Relative Dispersion:All --Interpolated');
            subplot(5,2,6);
            imagesc(FDmapPS_interp(:,:,a));colormap('jet'); caxis([1 1.5])
            colorbar; axis('image'); axis off;
            title('FD Map: Power Spectrum --Interpolated');
            subplot(5,2,7);
            imagesc(RobustFDStats(:,:,a,1));colormap('jet'); caxis([0 0.1]);
            colorbar; axis('image'); axis off;
            title('P Value of fit - RD Fast');
            subplot(5,2,8);
            imagesc(RobustFDStats(:,:,a,2));colormap('jet'); caxis([0 0.1]);
            colorbar; axis('image'); axis off;
            title('P Value of fit - RD Slow');
            subplot(5,2,9);
            imagesc(RobustFDStats(:,:,a,3));colormap('jet'); caxis([0 0.1]);
            colorbar; axis('image'); axis off;
            title('P Value of fit - RD All');
            subplot(5,2,10);
            imagesc(RobustFDStats(:,:,a,4));colormap('jet'); caxis([0 0.1]);
            colorbar; axis('image'); axis off;
            title('P Value of fit - FD Power');
            % save images as colour postscript images at 1200 dpi
            eval(['print -dpsc -r1200 ' num2str(FDdir) '_FD_Slice' num2str(a+Za-1) ' ;']); %AMW
            close %AMW
        end
        
    elseif doComponents=='n'
        
        for a=1:ZZ;  % % MoAl_Jan2011 %AMW
            figure(a); %AMW
            subplot(3,2,1);
            imagesc(bold_mean(:,:,a)); axis('image'); colormap('gray'); axis off;
            title('BOLD MEAN');
            subplot(3,2,2);
            imagesc(bold_std(:,:,a)); axis('image'); colormap('gray'); axis off;
            title('BOLD STDev');
            subplot(3,2,3);
            imagesc(FDmapRDall_interp(:,:,a));colormap('jet');
            colorbar; axis('image'); axis off;
            title('FD Map: Relative Dispersion');
            subplot(3,2,4);
            imagesc(FDmapPS_interp(:,:,a));colormap('jet'); caxis([1 1.5])
            colorbar; axis('image'); axis off;
            title('FD Map: Power Spectrum');
            subplot(3,2,5);
            imagesc(RobustFDStats(:,:,a,3));colormap('jet'); caxis([0 0.1]);
            colorbar; axis('image'); axis off;
            title('P Value of fit - RD All');
            subplot(3,2,6);
            imagesc(RobustFDStats(:,:,a,4));colormap('jet'); caxis([0 0.1]);
            colorbar; axis('image'); axis off;
            title('P Value of fit - FD Power');
            eval(['print -dpsc -r1200 ' num2str(FDdir) '_FD_Slice' num2str(a+Za-1) ' ;']); %AMW
            close %AMW
        end
    end
end;

disp('Finished saving FD maps');
toc
tic

% ----------------- OPTIONAL Save Figures [END]-----------

% ----------------- OPTIONAL ROI analysis [START]-----------

% % CREATING ROI MATRIX from MRS cooridinates % needs MRSx and MRSy to be
% entered before run this program
if automation=='n'              % MAW
    corin=0;
    while (corin==0) %AMW
        disp(' ')
        ROIyes=input('Will this analysis include ROIs yes/no (y/n) ? ','s');
        if ROIyes=='y' || ROIyes=='n'
            corin=1;
        else
            display('try again silly')
        end
    end
end

if ROIyes=='y',
    
    if automation=='n'              % MAW
    
        corin=0;
        while (corin==0) %AMW
            ROItype=input('Would you like to Draw an ROI, import a Mask or input Coordinates (d/m/c)? ','s');
            if ROItype=='c' || ROItype=='m' || ROItype=='d'
                corin=1;
            else
                display('Come on now...try again')
            end
        end
        
    end
    
    if ROItype=='c' %AMW
        if automation=='n'              % MAW
            ROIname=input('What would you like to name the ROI? ','s'); %AMW
            rn=input('What is the Right-most coord? (note: R is negative) '); %AMW
            ln=input('What is the Left-most coord? (note: R is negative) '); %AMW
            an=input('What is the Anterior-most coord? (note: A is negative) '); %AMW
            pn=input('What is the Posterior-most coord? (note: A is negative) '); %AMW
            in=input('What is the Inferior-most coord? (note: I is negative) '); %AMW
            sn=input('What is the Superior-most coord? (note: I is negative) '); %AMW
        end
        xrnd=round((ln-rn)/Info.DELTA(1)); %AMW
        yrnd=round((pn-an)/Info.DELTA(2)); %AMW
        zrnd=round(((sn-in)/Info.DELTA(3))+1); % AMW: add one, because this is used to decide how many slices
        MRSx=(rn+ln)/2; %AMW
        MRSy=(an+pn)/2; %AMW
        roioriginx=uint8(((Info.ORIGIN(1)*(-1)+MRSx*(-1))/Info.DELTA(1))-(xrnd/2)); %AMW
        roioriginy=uint8(((Info.ORIGIN(2)*(-1)+MRSy*(-1))/Info.DELTA(2))-(yrnd/2)); %AMW
        block1=zeros(roioriginx,roioriginy,'uint8');
        block2=ones(xrnd,yrnd,'uint8'); %AMW
        block3=zeros((Info.DATASET_DIMENSIONS(1)-roioriginx-xrnd),(Info.DATASET_DIMENSIONS(2)-roioriginy-yrnd),'uint8'); %AMW
        % SLnum=uint8((Info.ORIGIN(3)-sn)/(Info.DELTA(3)*(-1))):uint8(((Info.ORIGIN(3)-sn)/(Info.DELTA(3)*(-1)))+zrnd); %AMW
        lowsln=uint8((Info.ORIGIN(3)-sn)/(Info.DELTA(3)*(-1)));
        highsln=uint8(((Info.ORIGIN(3)-sn)/(Info.DELTA(3)*(-1)))+zrnd);
        
        corin=0;
        while (corin==0)
            %disp([''])
            SLnum=input(['Pick a slice between: ' num2str(lowsln) ' and ' num2str(highsln) ' ']);
            if SLnum <= highsln && SLnum >= lowsln
                corin=1;
            else
                disp('Your no good at this, try again')
                disp(' ')
            end
        end
        BW=blkdiag(block1,block2,block3);
        BW=double(BW);
        
    elseif ROItype=='d'
        
        % display mean BOLD map for doing ROI selection eventually this
        % will need replaced with high res anatomic anatomic data needs
        % registered to BOLD
        SLnum=input('What slice would you like to draw on? '); %AMW
        % SLnum=(SLnum-realZZ)*(-1); %slice number, as compared to AFNI, is backwards...slice 0 in AFNI is slice 30 in this program... this line fixes that error
        figure(1);
        imagesc(bold_mean(:,:,SLnum)); brighten(0.3); axis off; axis image; colormap('gray');
        disp('RIGHT mouse clicks for ROI, then LEFT mouse button for end');
        
        for dji=1:100,
            [x,y,BW,xi,yi]=roipoly;
            xi=floor(xi);
            yi=floor(yi);
            hold on;
            plot(xi,yi,'r');
            likeROI=input('What do you think ... Is this ROI a keeper ? ','s');
            if likeROI == 'y'
                ROIname=input('Give a name to this ROI ','s');
                % save image showing ROI location
                eval(['print -f1 -dpsc -r1200 ' num2str(FDdir) 'roi_' ROIname ';']);
                break
            else
                delete(1); figure(1);
                imagesc(bold_mean(:,:,SLnum)); brighten(0.3); axis off; axis image; colormap('gray');
                clear x y BW xi yi;
            end;
        end;
        
        BW=double(BW);  % Convert from uint8 to double precision for calculations
        newROIsize = nnz(BW);  % number of nonzero elements;
        
    elseif ROItype=='m'     % % MAW_Jan2011
        
        % % MAW_Jan2011: mask should be an AFNI file that is single sliced
        % and matched the resolution of the BOLD data % this can be done
        % using 3dresample eg: 3dresample -master FD_6447+orig -prefix
        % mask_resample_6447 -inset mask_6447+orig
        
        disp(' ');
        
        if automation=='n'
            % maskfile=input('Mask file (e.g. mask+orig): ','s');
            
            
            disp(' ')
            disp('TYPE ANY KEY to continue to select the mask file you would like to open (BRIK file)')
            
            pause
            maskfile=uigetfile('*.BRIK', 'Select the BRIK file for your mask'); % % MAW
        end
        
        maskfile = sprintf('%s', maskfile);   % % convert to string
        [mask_err,MASK,mask_Info,mask_ErrMessage]=BrikLoad(maskfile);
        % % MAW_Jan2011 imports mask as matrix with ROIs (value is 1, 2,
        % 3...) which can be separated into three matrixes if needed
        % begin with assumption that there is only one ROI (value would be
        % = 1)
        MASK=permute(MASK,[2,1]);
        
        disp(' ')
        disp(['Your data has ' num2str(max(max(MASK))) ' ROI(s).']); % MAW_Feb2011
        
        if max(max(MASK))>1
            
            if automation=='n'
                corin=0;
                while (corin==0) %MAW
                    roi_number=input('Which ROI would you like to analyse (1,2 etc)? '); % MAW_Feb2011
                    if roi_number <= num2str(max(max(MASK)))
                        corin=1;
                    else
                        display ('Try again Bonzo the clown')
                    end
                end
                
            end
            
            roi_mask=zeros(64,64);
            for a=1:64
                for b=1:64
                    if MASK(a,b)==roi_number
                        roi_mask(a,b)=1;
                    end
                end
            end
            
        end
        
        if automation=='n'
            
            disp ('Which slice did your mask come from');
            SLnum=input('Pick a slice: ');
                
            ROIname=input('What would you like to name the ROI? ','s'); %MAW
            
        end
        
        roi_mask=MASK;
        
        BW=double(roi_mask); % % assuming mask is single slice
        
    end
    
    BW_interp=imresize(BW,'scale', 8, 'method', 'nearest'); % % MAW_Feb2010 used for saving ROI data
    
    % Multiply all non-ROI elements with zero to create a mask;
    
    if doComponents=='y'
        ROIFDmapRDslow = BW.*FDmapRDslow(:,:,SLnum);
        ROIFDmapRDfast = BW.*FDmapRDfast(:,:,SLnum);
        ROI_fast_interp = BW_interp.*FDmapRDfast_interp(:,:,SLnum);    % % MAW_Feb2010
        ROI_slow_interp = BW_interp.*FDmapRDslow_interp(:,:,SLnum);    % % MAW_Feb2010
    end
    
    ROIbold_mean = BW.*bold_mean(:,:,SLnum);
    ROIbold_std = BW.*bold_std(:,:,SLnum);
    ROIFDmapPS = BW.*FDmapPS(:,:,SLnum);
    ROIFDmapRDall = BW.*FDmapRDall(:,:,SLnum);      % % MAW_Feb2010
    ROI_all_interp = BW_interp.*FDmapRDall_interp(:,:,SLnum);    % % MAW_Feb2010
    ROI_PS_interp = BW_interp.*FDmapPS_interp(:,:,SLnum);    % % MAW_Feb2010
    
    % locate all non-zero values inside ROI and assign into new matrices
    
    if doComponents=='y'
        [qhy3,qhx3,FDRDslow]=find(ROIFDmapRDslow);
        [qhy4,qhx4,FDRDfast]=find(ROIFDmapRDfast);
        [qhy7,qhx7,ROI_fast_interp]=find(ROI_fast_interp);        % % MAW_Feb2010
        [qhy8,qhx8,ROI_slow_interp]=find(ROI_slow_interp);        % % MAW_Feb2010
    end
    
    [qhy1,qhx1,MN]=find(ROIbold_mean);
    [qhy2,qhx2,stdev]=find(ROIbold_std);
    [qhy5,qhx5,FDps]=find(ROIFDmapPS);
    [qhy6,qhx6,FDRDall]=find(ROIFDmapRDall);        % % MAW_Feb2010
    [qhy9,qhx9,ROI_all_interp]=find(ROI_all_interp);        % % MAW_Feb2010
    [qhy10,qhx10,ROI_PS_interp]=find(ROI_PS_interp);        % % MAW_Feb2010
    
    % ROI SUMMARY STATISTICS ---------
    %   Testing for normal distribution
    %      H = 0 => Do not reject the null hypothesis at significance level
    %      ALPHA. H = 1 => Reject the null hypothesis at significance level
    %      ALPHA.
    %     The Bera-Jarque hypotheses are: Null Hypothesis:        X is
    %     normal with unspecified mean and variance. Alternative
    %     Hypothesis: X is not normally distributed.
    if doComponents=='y'
        mean_FDRDslow = mean(FDRDslow); median_FDRDslow = median(FDRDslow); std_FDRDslow = std(FDRDslow); jbtest_FDRDslow = jbtest(FDRDslow,0.01);
        mean_FDRDfast = mean(FDRDfast); median_FDRDfast = median(FDRDfast); std_FDRDfast = std(FDRDfast); jbtest_FDRDfast = jbtest(FDRDfast,0.01);
        sig_fit_fast = sum(sum(sum(RobustFDStats(:,:,:,1)<=0.05)));
        sig_fit_slow = sum(sum(sum(RobustFDStats(:,:,:,2)<=0.05)));
    end
    
    mean_BOLD = mean(MN); median_BOLD = median(MN); std_BOLD = std(MN); jbtest_BOLD = jbtest(MN,0.01);
    mean_stdev = mean(stdev); median_stdev = median(stdev); std_stdev = std(stdev); jbtest_stdev = jbtest(stdev,0.01);
    mean_FDRDall = mean(FDRDall); median_FDRDall = median(FDRDall); std_FDRDall = std(FDRDall); jbtest_FDRDall = jbtest(FDRDall,0.01);
    mean_FDps = mean(FDps); median_FDps = median(FDps); std_FDps = std(FDps); jbtest_FDps = jbtest(FDps,0.01);
    if ROIyes=='y'
        RobustFDStats(:,:,SLnum,3)= BW.*RobustFDStats(:,:,SLnum,3);
        RobustFDStats(:,:,SLnum,4)= BW.*RobustFDStats(:,:,SLnum,4);
        ROI_size=sum(sum(BW));
        sig_fit_all = ROI_size-(sum(sum(sum(RobustFDStats(:,:,SLnum,3)>0.05))));
        sig_fit_PS = ROI_size-(sum(sum(sum(RobustFDStats(:,:,SLnum,4)>0.05))));
    elseif ROIyes=='n'
        sig_fit_all = sum(sum(sum(RobustFDStats(:,:,:,3)<=0.05)));
        sig_fit_PS = sum(sum(sum(RobustFDStats(:,:,:,4)<=0.05)));
    end
        
    
    % % <= will give a value of 1 for each value where it is true. % sum
    % will add values for columns then sum does the same for rows %
    % effectively we get the number of voxels with sig fits % usually 4096
    % voxels total in ROI (64X64)
    % create a 2D variable called statsdata that contains all the stats
    % calculated above;
    
    
    
    if doComponents=='y'
        statsdata=[mean_BOLD,median_BOLD,std_BOLD,jbtest_BOLD;
            mean_stdev,median_stdev,std_stdev,jbtest_stdev;
            mean_FDRDslow,median_FDRDslow,std_FDRDslow,jbtest_FDRDslow;
            mean_FDRDfast,median_FDRDfast,std_FDRDfast,jbtest_FDRDfast;
            mean_FDRDall,median_FDRDall,std_FDRDall,jbtest_FDRDall;
            mean_FDps,median_FDps,std_FDps,jbtest_FDps;
            sig_fit_fast,sig_fit_slow,sig_fit_all,sig_fit_PS];
    else
        statsdata=[mean_BOLD,median_BOLD,std_BOLD,jbtest_BOLD;
            mean_stdev,median_stdev,std_stdev,jbtest_stdev;
            mean_FDRDall,median_FDRDall,std_FDRDall,jbtest_FDRDall;
            mean_FDps,median_FDps,std_FDps,jbtest_FDps;
            ROI_size,0,sig_fit_all,sig_fit_PS];
    end
    
    % save the stats as a space delimited 8-digit ASCII file to load into
    % excel, etc.
    ROIname=char(ROIname);
    eval(['save ' num2str(FDdir) '_' ROIname '_slice_' SLnum '_roistats.txt' ' statsdata -ASCII;']);
    clear mean_BOLD median_BOLD std_BOLD jbtest_BOLD;
    clear mean_stdev median_stdev std_stdev jbtest_stdev;
    clear mean_FDRDslow median_FDRDslow std_FDRDslow jbtest_FDRDslow;
    clear mean_FDRDfast median_FDRDfast std_FDRDfast jbtest_FDRDfast;
    clear mean_FDRDall median_FDRDall std_FDRDall jbtest_FDRDall;
    clear mean_FDps median_FDps std_FDps jbtest_FDps;
    clear sig_fit_fast sig_fit_slow sig_fit_all sig_fit_PS;
    
    % END of ROI selection loop----
end

disp('Finished ROI stats');
toc

% So-- save data to a matlab .mat file for later
% AMW, radically simplified this... Further commented out because, I really
% don't see the need for it.

% if ROIyes=='y'
%     if doComponents=='y' %AMW add filename for the kind of ROIstats...
%         if doSlice=='a' && bins==powertwo-2
%             eval(['save ' num2str(FDdir) '_ROIstats_all_Components statsdata ROI_all_interp;']);
%         else
%             eval(['save ' num2str(FDdir) '_ROIstats_Slices_' num2str(Za) '-' num2str(Zb) '_Bins_' num2str(bins) '_Components statsdata ROI_all_interp;']);
%         end
%     else
%         if doSlice=='a' && bins==powertwo-2
%             eval(['save ' num2str(FDdir) '_ROIstats_all statsdata ROI_all_interp;']);
%         else
%             eval(['save ' num2str(FDdir) '_ROIstats_Slices_' num2str(Za) '-' num2str(Zb) '_Bins_' num2str(bins) ' statsdata ROI_all_interp;']);
%         end
%     end
% end

% ----------------- OPTIONAL ROI analysis [END]-----------

%=========================================================

disp('Okay.... ready to do the next one!');
% clear  % % removed this for when this program is run in a loop
