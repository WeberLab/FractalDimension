function result=linreg(x,y);
% function result=linreg(x,y);
% Calculates the slope and the correlation coefficent of fitted line
% used by linear regression
% y=ax+b;
% result(1): slope (a);
% result(2): coefficient (r);
%
% WRITTEN BY:
% Dr. Peter Herman
% Dr. Andras Eke
% Fractal Physiology Lab
% Experimental Research Department
% II. Institute of Physiology
% Budapest, Ulloi ut 78-A
% Hungary 1082
%
% RIGHTS AND ACKNOWLEDGEMENTS:
% This software was developed with support from OTKA Grants I/3 2040, T 016953, 
% and NIH grant TW00442. Please cite these grants in any publication for which 
% this software or any derivates of it is used and send one reprint to Dr. Eke
% at the address given above.

mikro=mean(x);
nu=mean(y);
Qxy=sum((x-mikro).*(y-nu));
Qx=sum((x-mikro).^2);
Qy=sum((y-nu).^2);
result(1)=Qxy/Qx;
result(2)=Qxy/(Qx*Qy)^0.5;