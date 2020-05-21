function result=bridge(ts);
% function result=bridge(ts);
% Bridge deterended method: input ts, output BD ts
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

n=length(ts);
points=n-1:-1:0;
line=((ts(1)-ts(n))*points/(n-1))+ts(n);
line=line'; 
result=ts-line;
