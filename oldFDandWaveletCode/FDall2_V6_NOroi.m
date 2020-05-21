% M. Noseworthy, PhD, & G.M. Wardlaw, MSc. v10.9 
% IMPORTANT- AFNI Matlab toolbox is required for THIS code.
%            get this from: http://afni.nimh.nih.gov/afni/matlab
%            This also includes info on how HEAD/BRIK load/save works

%  HERE's what this program assumes you have already done:
%  1.  The time series data has been renamed to I.#.# scheme using AFNIrename
%  2.  You have used AFNI "to3d" to save BOLD data as HEAD and BRIK files
% e.g. to3d -epan -time:tz 2400 3 250 altplus -save_outliers bad.txt Ser3/I.0*.*
%  3.  Make sure you run this .m file in the directory with all your AFNI data!!!

% NOTE this uses MATLAB AFNI toolbox functions 'BrikLoad' and 'WriteBrik';

% ========================================================================
% MDN Sept 9/09  Used permute command to rotate matrix dimensions (thx Arv)
% MDN Sept 9/09  Added functionality to work with any sie BOLD matrices

% ========================================================================

% % ================================================================
% % MAW 10 NOV 2009 modified so that defaults are set:
% % Must enter 'examID=XXXX' 'MRSx=XX.X' 'MRSy=XX.X' before running
% % Defaults: SLnum = 1, accept ROI (matrix made from MRS voxel)
% % FDdo (previous FD?) set to 'n' which means it has to DO FD calc
% %
% % This was done to allow batch running of this program
% %
% % ================================================================

% % MAW 14 Dec 209 V4 changes MRSy to + from - since the input Y 
% % coordinates now have the correct sign.

% % MAW_Feb2010
% % V5 now takes 2048 data points from the MIDDLE of the data set
% % Also, added stats for RD and PS fits, added a new matrix of these stats
% % Added a new RD fit, with all data (not assume two trends)
% % changed the way the figures are produced
% % changed the stats to include the number of significant fits

% % MAW_June2010
% % modified so that MRSx and MRSy do not need to be entered, and exam ID
% % is asked for. This way, no variables have to be entered in advance. 
% % Also, ROI analysis has been turned of so that all we get is FD maps.
% % This can be cleaned up later to completely remove the ROI loop.

% % Note: comments added by MAW usually begin with '% %'

tic
beep off;  % I hate this fucking noise when I have my headphones on!

disp('This program is automated multislice FD analysis saving as Matlab HEAD/BRIK files');
disp(' ');
disp('This program has 3 parts:');
disp('   1) Creates .mat datafiles of various FD parameters');
disp('   2) Save parametric FD maps as AFNI HEAD/BRIK fbuc files'); %%%%NOT DONE
%%%%% #2 NOT DONE- see line 223!
disp('   3) Saves high res postscript files of FD parametric maps');
disp('   4) Performs ROI analysis');
disp('This program can work on previously saved .mat file of FD data');

% % clear all;
close all;
warning off;

ls

MRSx=20;    % % default values if no ROI is being used
MRSy=20;    % % default values if no ROI is being used


FDdir=input('What is the 4 digit exam number?','s');

examID=FDdir;

eval(['cd FD_' num2str(FDdir) ';']);
ls
disp(' ');
disp('There may already be some FD maps in this directory');
% % doFD=input('Is there a .mat FDmap file already saved (y/n) ','s');
doFD='n';

if doFD=='n',
    
% Input variables
% % Bbold=input('BOLD DATA (e.g. data+orig): ','s');  
Bbold = ['bold_',num2str(examID),'_FD+orig'];

% % TR=input('Repitition Time (TR) for BOLD data = ');
TR=0.25;

% Load BOLD data BRIK using AFNI plugin 'brikload'
% This will create a 4D variable, BOLD, with xres x yres x #slices x time
% NOTE- will need to keep array 'Info' as this is the header info and can be used
% for writing BRIK files later

brainbold = sprintf('%s', Bbold);   % % convert to string
[err,BOLD,Info,ErrMessage]=brikload(brainbold);

% ########## TRANSPOSE DATA FOR CORRECT DISPLAY ORIENTATION ############################
% This rotates images 90 so they align up-down
% so BOLD data in correct orientation: (LR,AP,SI,time)
%newBOLD=zeros(XX,YY,ZZ,2400);
BOLD=permute(BOLD,[2,1,3,4]);   % % Rearrange dimensions of array

XX=size(BOLD,1);
YY=size(BOLD,2);
ZZ=size(BOLD,3);


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
% ------------ Perform FD calculations on all data [START]------------------
% Here we are fitting 2 components- the fast and slow FDs.

%Assigning zero filled matrices for FD RD fast, slow and FD PS
FDmapRDall=zeros (XX,YY,ZZ);  % % MAW_Feb2010
FDmapRDfast=zeros(XX,YY,ZZ);
FDmapRDslow=zeros(XX,YY,ZZ);
FDmapPS=zeros(XX,YY,ZZ);
RobustFDStats=zeros (XX,YY,4); % % MAW_Feb2010 - 3rd dimension for Stats of RDfast, RDslow, RDall and FD_power (PS)

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
    
	fractaldataRD=zeros(10,2);
	nosigcount=0;
    	% note- assumes BOLD is done with discarded acquisitions (re: steadystate).
	% Ztemp=BOLDtemp(r,c,1:2048);  % takes the first 2048 time pts out of 2400 for slice selection
    Ztemp=BOLDtemp(r,c,177:2224);     % % MAW_Feb2010 takes the middle 2048 time pts out of 2400 for slice selection
    
	Z=(squeeze(Ztemp))'; %extract and transpose vector for pixel location 1,1.
	
	sumcheck=sum(Z,2);	

		if sumcheck == 0, 
			% *** NO MR SIGNAL, ASSUME FD is NOISE ==1.5 *****
			FDmapRDfast(r,c,SL) = 1.5;
			FDmapRDslow(r,c,SL) = 1.5;
            FDmapRDall (r,c,SL) = 1.5;  % % MAW_Feb2010
            RobustFDStats(r,c,:) = 0;   % % MAW_Feb2010
  		     	clear sumcheck;
			nosigcount=nosigcount+1;
		else
			for d=0:9, % i.e. 2^9=512       % % All data points are merged into sucessively smaller number of data points by averaging adjacent data points (i.e. they are binned)
				m=2^d;                      % % max number of points averaged is 512 which leaves us with 4 bins
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
                        
			FD_RDfast = fractaldataRD(1:4,:); % 1st 4 points used in FDfast fit     % % more data points, smaller bins ('m' is low)
			FD_RDslow = fractaldataRD(5:9,:); % 1ast 5 points used in FDslow fit    % % fewer data points because large bins have been averaged
            FD_RDall = fractaldataRD; % % MAW_Feb2010 - all RD time points
            

            % Perform fast and slow (and all) component assessment
            
            xmfast=log10(FD_RDfast(:,1));   % % first column is 'm' which was number in each bin
			ymfast=log10(FD_RDfast(:,2));   % % second column was the RD
			xmslow=log10(FD_RDslow(:,1));
			ymslow=log10(FD_RDslow(:,2));
            xmall=log10(FD_RDall(:,1)); % % MAW_Feb2010
            ymall=log10(FD_RDall(:,2)); % % MAW_Feb2010
			
            % Check for NaN values using function 'isnan'
            if sum(isnan(xmfast))==0 && sum(isnan(ymfast))==0 && sum(isnan(xmslow))==0 && sum(isnan(ymslow))==0 && sum(isnan(xmall))==0 && sum(isnan(ymall))==0,
                [fitsfast, stats_fast]=robustfit(xmfast,ymfast);      % % MAW_Feb2010 added stats meassure to each fit
                FDmapRDfast(r,c,SL)=1-(fitsfast(2,1));
                RobustFDStats(r,c,1) = stats_fast.p(2,1);            % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 1 for RDfast)
                                                            % % Note that stats.p gives P value for the y-intercept and for the slope (2X1 matrix)
            	[fitsslow, stats_slow]=robustfit(xmslow,ymslow);
                FDmapRDslow(r,c,SL)=1-(fitsslow(2,1));
                RobustFDStats(r,c,2) = stats_slow.p(2,1);            % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 2 for RDslow)
                
                [fitsall, stats_all]=robustfit(xmall,ymall);         % % MAW_Feb2010 - robustfit is 2X1 (y-intercept, slope)
                FDmapRDall(r,c,SL)=1-(fitsall(2,1));    % % MAW_Feb2010 - FD(RD) = 1-slope 
                RobustFDStats(r,c,3) = stats_all.p(2,1);             % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 3 for RDall)
                
            else
                FDmapRDfast(r,c,SL)=1.5;
                FDmapRDslow(r,c,SL)=1.5;
                FDmapRDall(r,c,SL)=1.5;     % % MAW_Feb2010
                RobustFDStats(r,c,:) = 0;   % % MAW_Feb2010
            end;
            
	


    
% NOW finished making: FDmapRDfast, FDmapRDslow for slice=SL
% % and FDmapRDall


% ####  POWER SPEC METHOD  -----------------------------------
		fractaldataPS=zeros(10,2);
		sumcheck=sum(Z,2);
		
		if sumcheck == 0,
			% *** NO MR SIGNAL, ASSUME==FD is NOISE ==1.5 *****
			FDmapPS(r,c,SL) = 1.5;
            RobustFDStats(r,c,:) = 0;   % % MAW_Feb2010
  		     	clear sumcheck;
			nosigcount=nosigcount+1;
		else
			Z=detrend(Z);           % % "removes the linear trend from a matrix, usually for FFT processing"
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
			D=2-hurst;	  %by selecting D=2- not 1-, scales similarly with RD.
				          %usually calculated as 2- from article.
			FDmapPS(r,c,SL)=D;
            RobustFDStats(r,c,4) = statsPS.p(2,1);   % % MAW_Feb2010 - add P value of fit to a matrix ("slice" 4 for PS)
		end;
	            
	c=c+1; %column counter-----
	end;  

r=r+1; %row counter-----
end;  

clear Ztemp;
clear Z;
clear sumcheck;

%Graeme: Clear BOLDtemp, and re-assign for slice SL=1+SL at top of loop (line 110)
SL=SL+1;
clear BOLDtemp;

end;

disp('Starting image interpolation');
toc
tic

% DO These interpolations After Loop!!
FDmapRDfast_interp=imresize(FDmapRDfast,[interpsize interpsize],'bicubic');
FDmapRDslow_interp=imresize(FDmapRDslow,[interpsize interpsize],'bicubic');
FDmapRDall_interp=imresize(FDmapRDall,[interpsize interpsize],'bicubic');
FDmapPS_interp=imresize(FDmapPS,[interpsize interpsize],'bicubic');


% ------------ Perform FD calculations [END]-----------


% So-- save data to a matlab .mat file for later
eval(['save ' num2str(FDdir) '_FDmaps FDdir Info XX YY ZZ bold_mean bold_std FDmapRDslow FDmapRDfast FDmapRDall FDmapPS RobustFDStats FDmapRDslow_interp FDmapRDfast_interp FDmapRDall_interp FDmapPS_interp;']);
% save temporary file sss.  This has only the directory name
%save sss FDdir



% ----------------- Create AFNI HEAD/BRIK files [START]-----------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is where a section to write AFNI HEAD/BRIK files needs to be added
% To do this you will need to make use of AFNI matlab function 'Writebrik'

% 1)  You will need: 'Info' made by AFNI plugin 'brikload' (top of this pgm) 
% 2)  The creation of BRIK/HEAD files of the following wlll be required: 
%     bold_mean bold_std FDmapRDslow FDmapRDfast FDmapPS 
% 3)  These are all 3D variables AP,LR,SI- or 64x64 x 3 slices
% 4)  There is no need to make BRIK/HEAD for interpolated versions of these 
% as this can be done using the AFNI GUI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------- Create AFNI HEAD/BRIK files [END]-----------



else
    eval(['load ' num2str(FDdir) '_FDmaps ;']);    % % from IF statement on line 65
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parts below are optional and prompt the user for decision as to whether
% they are to be done or not.

% ----------------- OPTIONAL Save Figures [START]-----------
% % gofig=input('Do you want to make high res ps image files ? ','s');
gofig='y';

if gofig == 'y',
    
    
% Display Mean of each of the 3 slices and pick one to analyze -----------
for n=1:ZZ;
   figure(1);
   subplot(1,3,n); imagesc(bold_mean(:,:,n)); 
   axis('image'); colormap('gray'); axis off;
   title(['Mean BOLD: ',num2str(n)]);
   hold on; 
end;
    
% % SLnum=input('Pick a slice:-- ');
SLnum=1;

delete(1);


% Create parametric maps from each slice-------------
figure(1);
    subplot(1,2,1);         % % Sub divide figure in a 1 X 2 'sub plots' and place this one in box 1 
        imagesc(bold_mean(:,:,SLnum)); axis('image'); colormap('gray'); axis off;  
        title('BOLD MEAN');
    subplot(1,2,2);
        imagesc(bold_std(:,:,SLnum)); axis('image'); colormap('gray'); axis off;  
        title('BOLD STDev');

figure(2);
    subplot(2,2,1);         % % Note, fills sub plots across row then down to next column
        imagesc(FDmapRDfast_interp(:,:,SLnum));colormap('jet');
        colorbar; axis('image'); axis off;
        title('FD Map: Relative Dispersion:FAST --Interpolated');
    subplot(2,2,2);
        imagesc(FDmapRDslow_interp(:,:,SLnum));colormap('jet');
        colorbar; axis('image'); axis off;
        title('FD Map: Relative Dispersion:SLOW --Interpolated');
    subplot(2,2,3);
        imagesc(FDmapRDall_interp(:,:,SLnum));colormap('jet');
        colorbar; axis('image'); axis off;
        title('FD Map: Relative Dispersion:All --Interpolated');
    subplot(2,2,4);
        imagesc(FDmapPS_interp(:,:,SLnum));colormap('jet'); caxis([1 1.5])
        colorbar; axis('image'); axis off;
        title('FD Map: Power Spectrum --Interpolated');

figure(3);
    subplot(2,2,1);
        imagesc(RobustFDStats(:,:,1));colormap('jet'); caxis([0 0.1]);
        colorbar; axis('image'); axis off;
        title('P Value of fit - RD Fast');
    subplot(2,2,2);
        imagesc(RobustFDStats(:,:,2));colormap('jet'); caxis([0 0.1]);
        colorbar; axis('image'); axis off;
        title('P Value of fit - RD Slow');
    subplot(2,2,3);
        imagesc(RobustFDStats(:,:,3));colormap('jet'); caxis([0 0.1]);
        colorbar; axis('image'); axis off;
        title('P Value of fit - RD All');
    subplot(2,2,4);
        imagesc(RobustFDStats(:,:,4));colormap('jet'); caxis([0 0.1]);
        colorbar; axis('image'); axis off;
        title('P Value of fit - FD Power');
        

        
disp('Maximize the sizes of the images!!');  % enhances images resolution in final image files;
disp('TYPE THE ANY KEY to continue');
% % pause

% save images as colour postscript images at 1200 dpi
eval(['print -f1 -dpsc -r1200 ' num2str(FDdir) '_Mean_StDev_V2 ;']);
eval(['print -f2 -dpsc -r1200 ' num2str(FDdir) '_FD_Maps_V2 ;']);
eval(['print -f3 -dpsc -r1200 ' num2str(FDdir) '_P_Value_of_Fit_V2 ;']);

end;

disp('Finished saving FD maps');
toc
tic

% ----------------- OPTIONAL Save Figures [END]-----------

% % CREATING ROI MATRIX from MRS cooridinates 
% % needs MRSx and MRSy to be entered before run this program

roioriginx=uint8(((Info.ORIGIN(1)*(-1)+MRSx)/3.75)-3);

roioriginy=uint8(((Info.ORIGIN(2)*(-1)+MRSy)/3.75)-3);

block1=zeros(roioriginy,roioriginx,'uint8');
block2=ones(6,6,'uint8');
block3=zeros((64-roioriginy-6),(64-roioriginx-6),'uint8');
MRSroi=blkdiag(block1,block2,block3);
        
        
% ----------------- OPTIONAL ROI analysis [START]-----------

% % ROIyes=input('Will this analysis include ROIs (y/n) ? ','s');
ROIyes='n'; % % default value if no ROI is being used

% % if ROIyes=='y',
% %     
% % % Display Mean of each of the 3 slices and pick one to analyze -----------
% % for n=1:ZZ;
% %    figure(1);
% %    subplot(1,3,n); imagesc(bold_mean(:,:,n)); 
% %    axis('image'); colormap('gray'); axis off;
% %    title(['Mean BOLD: ',num2str(n)]);
% %    hold on;  
% % end;
    
% % SLnum=input('Pick a slice:-- ');
SLnum=1;

% % delete(1);

% So at this point we have 5 maps to get ROI stats from:
% cropbold_mean, cropbold_stdev, FDmapRDslow, FDmapRDfast, and FDmapPS
% % ROInum=input('How many ROIs do you want to analyze ? ');
ROInum=1;

for nums=1:ROInum,
% %     
% %     % display mean BOLD map for doing ROI selection
% %     % eventually this will need replaced with high res anatomic
% %     % anatomic data needs registered to BOLD 
% %     figure(1);
% %     imagesc(bold_mean(:,:,SLnum)); brighten(0.3); axis off; axis image; colormap('gray');
% %     disp('RIGHT mouse clicks for ROI, then LEFT mouse button for end');
% % 
% %   for dji=1:100,
% %     [x,y,BW,xi,yi]=roipoly;
% %     xi=floor(xi);
% %     yi=floor(yi);
% %     hold on;
% %     plot(xi,yi,'r');
% %         likeROI=input('What do you think ... Is this ROI a keeper ? ','s');
% %             if likeROI == 'y'
% %                 ROIint=input('Give a name to this ROI ','s');
% %                 % save image showing ROI location
% %                 eval(['print -f1 -dpsc -r1200 ' num2str(FDdir) 'roi.' num2str(nums) ';']);
% %                 break
% %             else,
% %                 delete(1); figure(1);
% %                 imagesc(bold_mean(:,:,SLnum)); brighten(0.3); axis off; axis image; colormap('gray');
% %                 clear x y BW xi yi;
% %             end;
% %     end;

    BW=MRSroi;
    

    BW=double(BW);  % Convert from uint8 to double precision for calculations
    newROIsize = nnz(BW);  % number of nonzero elements;
    
    BW_interp=imresize(BW,'scale', 8, 'method', 'nearest'); % % MAW_Feb2010 used for svaing ROI data


    % Multiply all non-ROI elements with zero to create a mask;
    ROIbold_mean = BW.*bold_mean(:,:,SLnum);
    ROIbold_std = BW.*bold_std(:,:,SLnum);
    ROIFDmapRDslow = BW.*FDmapRDslow(:,:,SLnum);
    ROIFDmapRDfast = BW.*FDmapRDfast(:,:,SLnum);
    ROIFDmapPS = BW.*FDmapPS(:,:,SLnum);
    ROIFDmapRDall = BW.*FDmapRDall(:,:,SLnum);      % % MAW_Feb2010
    ROI_fast_interp = BW_interp.*FDmapRDfast_interp(:,:,SLnum);    % % MAW_Feb2010
    ROI_slow_interp = BW_interp.*FDmapRDslow_interp(:,:,SLnum);    % % MAW_Feb2010
    ROI_all_interp = BW_interp.*FDmapRDall_interp(:,:,SLnum);    % % MAW_Feb2010
    ROI_PS_interp = BW_interp.*FDmapPS_interp(:,:,SLnum);    % % MAW_Feb2010
    
    
    % locate all non-zero values inside ROI and assign into new matrices
    [qhy1,qhx1,MN]=find(ROIbold_mean);
    [qhy2,qhx2,stdev]=find(ROIbold_std);
    [qhy3,qhx3,FDRDslow]=find(ROIFDmapRDslow);
    [qhy4,qhx4,FDRDfast]=find(ROIFDmapRDfast);
    [qhy5,qhx5,FDps]=find(ROIFDmapPS);
    [qhy6,qhx6,FDRDall]=find(ROIFDmapRDall);        % % MAW_Feb2010
    [qhy7,qhx7,ROI_fast_interp]=find(ROI_fast_interp);        % % MAW_Feb2010
    [qhy8,qhx8,ROI_slow_interp]=find(ROI_slow_interp);        % % MAW_Feb2010
    [qhy9,qhx9,ROI_all_interp]=find(ROI_all_interp);        % % MAW_Feb2010
    [qhy10,qhx10,ROI_PS_interp]=find(ROI_PS_interp);        % % MAW_Feb2010
    


    % ROI SUMMARY STATISTICS ---------
    %   Testing for normal distribution
    %      H = 0 => Do not reject the null hypothesis at significance level ALPHA.
    %      H = 1 => Reject the null hypothesis at significance level ALPHA.
    %     The Bera-Jarque hypotheses are: 
    %     Null Hypothesis:        X is normal with unspecified mean and variance.
    %     Alternative Hypothesis: X is not normally distributed.
    mean_BOLD = mean(MN); median_BOLD = median(MN); std_BOLD = std(MN); jbtest_BOLD = jbtest(MN,0.01);
    mean_stdev = mean(stdev); median_stdev = median(stdev); std_stdev = std(stdev); jbtest_stdev = jbtest(stdev,0.01);
    mean_FDRDslow = mean(FDRDslow); median_FDRDslow = median(FDRDslow); std_FDRDslow = std(FDRDslow); jbtest_FDRDslow = jbtest(FDRDslow,0.01);
    mean_FDRDfast = mean(FDRDfast); median_FDRDfast = median(FDRDfast); std_FDRDfast = std(FDRDfast); jbtest_FDRDfast = jbtest(FDRDfast,0.01);
    mean_FDRDall = mean(FDRDall); median_FDRDall = median(FDRDall); std_FDRDall = std(FDRDall); jbtest_FDRDall = jbtest(FDRDall,0.01);
    mean_FDps = mean(FDps); median_FDps = median(FDps); std_FDps = std(FDps); jbtest_FDps = jbtest(FDps,0.01);
    sig_fit_fast = sum(sum(RobustFDStats(:,:,1)<=0.05)); sig_fit_slow = sum(sum(RobustFDStats(:,:,2)<=0.05)); sig_fit_all = sum(sum(RobustFDStats(:,:,3)<=0.05)); sig_fit_PS = sum(sum(RobustFDStats(:,:,4)<=0.05)); 
    % % <= will give a value of 1 for each value where it is true. 
    % % sum will add values for columns then sum does the same for rows
    % % effectively we get the number of voxels with sig fits
    % % usually 4096 voxels total in ROI (64X64)
    
    
% create a 2D variable called statsdata that contains all the stats calculated above;     
statsdata=[mean_BOLD,median_BOLD,std_BOLD,jbtest_BOLD;
    mean_stdev,median_stdev,std_stdev,jbtest_stdev;
    mean_FDRDslow,median_FDRDslow,std_FDRDslow,jbtest_FDRDslow;
    mean_FDRDfast,median_FDRDfast,std_FDRDfast,jbtest_FDRDfast;
    mean_FDRDall,median_FDRDall,std_FDRDall,jbtest_FDRDall;
    mean_FDps,median_FDps,std_FDps,jbtest_FDps;
    sig_fit_fast,sig_fit_slow,sig_fit_all,sig_fit_PS];
% save the stats as a space delimited 8-digit ASCII file to load into excel, etc. 
 eval(['save ' num2str(FDdir) '_roi_stats_V2.txt' ' statsdata -ASCII;']);
 
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
eval(['save ' num2str(FDdir) '_variables_V2 examID statsdata MRSx MRSy ROI_fast_interp ROI_slow_interp ROI_all_interp ROI_PS_interp FDdir Info bold_mean bold_std FDmapRDslow FDmapRDfast FDmapRDall FDmapPS RobustFDStats FDmapRDslow_interp FDmapRDfast_interp FDmapRDall_interp FDmapPS_interp;']);
% save temporary file sss.  This has only the directory name
%save sss FDdir


% % end;

        
% ----------------- OPTIONAL ROI analysis [END]-----------



%=========================================================

% go back a directory and change the name to indicate FD maps done;
% eval(['unix(''' 'mv ' num2str(FDdir) ' ' num2str(FDdir) 'FD'')' ';']);
cd ..
disp('Okay.... ready to do the next one!');
% clear  % % removed this for when this program is run in a loop



