function output=fractool(ts, outputfilename);
% function output=fractool(ts,outputfilename);
%
% *** FracTool ***
%  (Main Program)
%
% ABOUT THE NAME:
% FracTool is a mosaic name composed from Fractal Tool, a set of software tools
% designed and developed to aid time series analysis using the concept of
% statistical fractals.
%
% FILE:
% fractool.m (Version 1.0).
%
% LAST MODIFIED AT:
% 12/31 1999
%
% WRITTEN BY:
% Dr. Peter Herman (herman@elet2.sote.hu)
% Dr. Andras Eke (eke@elet2.sote.hu)
% Fractal Physiology Laboratory
% Experimental Research Department
% II. Institute of Physiology
% Semmelweis University, Budapest
% Budapest, Ulloi ut 78-A
% Hungary 1082
%
% RIGHTS AND CONDITIONS OF USE:
% This is freeware. Please place the following reference in any publication for which 
% this software or any derivates of it is used and send one reprint to Dr. Eke 
% at the address given above:
%
% A. Eke, P. Hermán, J. B. Bassingthwaighte, G. M. Raymond, D. B. Percival, 
% M. Cannon, I. Balla, and C. Ikrényi. Physiological time series: distinguishing
% fractal noises from motions. Pflügers Archiv European Journal of Physiology, 
% 439(4):403-415, 2000
%
% This software is permitted to be distributed only as a complete compressed package
% (see below). The README file must be included.
%
% ACKNOWLEDGEMENTS:
% FracTool was developed with support from OTKA Grants I/3 2040, T 016953, 
% and NIH grant TW00442. 
%
% CONTENTS OF FRACTOOL PACKAGE:
% README.txt
% Bridge.m
% Disper.m
% Bdswv.m
% Stdn.m
% Linreg.m
% Spec.m
% Window.m
% Fractool.m
%
% DESCRIPTION:
% FracTool is a software tool for fractal analysis of time series based on the 
% signal summation conversion method (SSC) described in the following 
% publication
%
% A. Eke, P. Hermán, J. B. Bassingthwaighte, G. M. Raymond, D. B. Percival, 
% M. Cannon, I. Balla, and C. Ikrényi. Physiological time series: distinguishing
% fractal noises from motions. Pflügers Archiv European Journal of Physiology, 
% 439(4):403-415, 2000
%
% FracTool determines 
% 1) the signal class of the time series according to the fGn/fBm model of 
%    random signals (fGn=fractional Gaussian noise, fBm= fractional Brownian 
%    motion), 
% 2) the Hurst coefficients for fGn or fBm time series.
%
% This program has been written in the MATLAB language. It runs in the MATLAB 
% environment on different platforms, like PC, Mac and Unix. It has been tested 
% under MATLAB v4.2 and v5.2. Further downward compatibility has not been 
% explored.
%
% INPUT FILE:
% The input file to the program is a text file containing the data of the time 
% series in a single column of the following structure:
%
% data1
% data2
% ...
% dataN
% 
% The length of the time series actually analyzed by the program is truncated 
% to the longest length of 2^n. We don't recommend using longer file, than 2^18
% datapoints.
%
% First, the MATLAB program is run. At the prompt, issuing the following command
% will load the input file to the input matrix to the program
%
% >> load filename.txt 
%
% The name "filename" is given to the matrix.
%
% EXECUTION OF PROGRAM:
% >> output=fractool(filename,outputfilename);
%
% OUTPUT FILE:
% The output from the program is matrix "output" that is automatically saved to 
% disk in the file "outputfilename" once the program has ended. If no "outputfilename" is
% given, the default name "fractool.txt" would be used. Its structure is as follows
%
% Length (2^n)   	 16.0000000000
% mean_of_data   	 -0.0000000000
% SD_of_data     	 13.5121865024
% class          	  2.0000000000
% beta           	  1.5501109967
% H_PSD          	  0.2750554983
% H_Disp         	 -1.0000000000
% H_bdSWV        	  0.3100622409
% H_fGn          	 -1.0000000000
% H_fBm          	  0.2925588696
%
% Symbols used in "fractool.txt":
%
% Length (2^n)   Length of the analyzed time series in power of 2
% mean_of_data   Mean of the time series
% SD_of_data     SD of the time series
% class          Signal type. 0 := fGn, 1 := fGn/fBm boundary, 2 := fBm, 
%					  3 := outside of fGn/fBm model 
% beta			  beta estimated by the low_PSD_w,e method
% H_PSD			  Hurst coefficient (H) estimated by the low_PSD_w,e method (-1 
%					  when not applicable)
% H_Disp         Hurst coefficient (H) estimated by the Dispersional method (Disp) 
%					  (-1 when not applicable)
% H_bdSWV        H estimated by the bridge detrended scaled windowed variance 
%                method (bdSWV)(-1 when not applicable)
% H_fGn          Hurst coefficient (H) calculated as the mean of H_Disp and H_PSD
%					  (-1 when not applicable)
% H_fBm          Hurst coefficient (H) calculated as the mean of H_bdSWV and H_PSD
%					  (-1 when not applicable)
%
% RUNTIME MESSAGES:
% The following messages are printed when the analysis is progressing properly.
%
% » output=fractool(filename);
%
% getready__time =
%
%    0.0600
%
%
% Elapsed__time =
%
%    8.5200
%
%
% WARNINGS:
% Error messages of MATLAB is used to warn user of potential problems.
%
% EXTERNAL FUNCTIONS CALLED:
%
% For estimating beta:
% Spec.m -> to estimate beta and H_PSD
% bridge.m -> endmatching
% window.m -> windowing
%
%     Input parameters:ts,p,avg
%     Output parameter: result
%
%     ts= time series
%     p= nth power of 2, where 2^n<= length of time series
% 		avg: if avg=0, then spectra all, if avg=1, then spectra avg;
% 		result(1): Hurst coefficient;
% 		result(2): correlation coefficient (r) from linear regression;
% 		result(3): beta from Fourier analysis;
% 		result(4): Hurst coefficient (high frequencies excluded);
% 		result(5): correlation coefficient (r) from linear regression (high frequencies excluded);
% 		result(6): beta from Fourier analysis (high frequencies excluded);
%
% For estimating Hurst coefficients:
% Bdswv.m -> to estimate H_bdSWV and HH
% Disper.m -> to estimate H_Disp
%
%     Input parameters:ts,p
%     Output parameter: result
%
%     ts= time series
%     p= nth power of 2, where 2^n<= length of time series
%     result(1)= Hurst coefficient
%     result(2)= correlation coefficient (r) from linear regression
%
%     Other functions called by bdswv.m and disper.m are:
%
%          bridge.m -> bridge detrending
%          linreg.m -> linear regression
%          stdn.m -> SD
%

H_PSD=-1;
H_Disp=-1;
H_bdSWV=-1;
H_fGn=-1;
H_fBm=-1;
t0=clock;
n=length(ts);
i=1;
while n>=2^i,
	i=i+1;
end;
p=i-1;
tsid=ts(1:2^p);
clear ts;
mean_of_data=mean(tsid);
SD_of_data=std(tsid);
n=length(tsid);
getready__time=etime(clock,t0)

% Calculates beta using lowPSDw,e

result=spec(bridge(window(tsid)),p,0);

if(result(6)<0.38) & (result(6)>-1),
   H_PSD=(result(6)+1)/2;
   temp=disper(tsid,p);
   H_Disp=temp(1);
   H_fGn=(H_PSD+H_Disp)/2;
   class=0;
elseif (result(6)>=0.38) & (result(6)<=1.04),
   ts=cumsum(tsid);
   temp=bdswv(ts,p);
   HH=temp(1);
   if HH<0.8,
      H_PSD=(result(6)+1)/2;
      temp=disper(tsid,p);
      H_Disp=temp(1);
      H_fGn=(H_PSD+H_Disp)/2;
      class=0;
   elseif (HH>0.8) & (HH<1),
      class=1;
  	else HH>1,
      H_PSD=(result(6)-1)/2;
      temp=bdswv(ts,p);
      H_bdSWV=temp(1);
      H_fBm=(H_PSD+H_bdSWV)/2;
      class=2;
   end;
elseif (result(6)>1.04) & (result(6)<3),
   H_PSD=(result(6)-1)/2;
   temp=bdswv(tsid,p);
   H_bdSWV=temp(1);
   H_fBm=(H_PSD+H_bdSWV)/2;
   class=2;
else
   class=3;
end; 

output(1)=p;         
output(2)=mean_of_data;   
output(3)=SD_of_data;     
output(4)=class;           
output(5)=result(6);			  
output(6)=H_PSD;			   
output(7)=H_Disp;          
output(8)=H_bdSWV;         
output(9)=H_fGn;          
output(10)=H_fBm;          

resultsstring=[
'Length (2^n)   ';
'mean_of_data   ';
'SD_of_data     ';
'class          ';
'beta           ';
'H_PSD          ';
'H_Disp         ';
'H_bdSWV        ';
'H_fGn          ';
'H_fBm          '];

 if nargin<2, outputfilename='FracTool.txt'; end;
fid=fopen(outputfilename,'w');
for i=1:10;
	fprintf(fid,'%s',resultsstring(i,1:15));
	fprintf(fid,'\t%14.10f\r\n',output(i));
end;
fclose(fid); 
Elapsed__time=etime(clock,t0)
