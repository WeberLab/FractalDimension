function result=linreg(x,y);
% function result=linreg(x,y);
% Calculates the slope and the correlation coefficent of a linear fit
% based on linear regression model
% y=ax+b;
% result(1): slope (a);
% result(2): coefficient (r);
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
% A. Eke, P. Herm�n, J. B. Bassingthwaighte, G. M. Raymond, D. B. Percival, 
% M. Cannon, I. Balla, and C. Ikr�nyi. Physiological time series: distinguishing
% fractal noises from motions. Pfl�gers Archiv European Journal of Physiology, 
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

mikro=mean(x);
nu=mean(y);
Qxy=sum((x-mikro).*(y-nu));
Qx=sum((x-mikro).^2);
Qy=sum((y-nu).^2);
result(1)=Qxy/Qx;
result(2)=Qxy/(Qx*Qy)^0.5;