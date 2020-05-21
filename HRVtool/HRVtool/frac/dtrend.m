function result=dtrend(x,y);
% function result=dtrend(x,y);
% x datokat a legkisebb modszer szerint detrendeli,

% y=a+bx fuggvenybol a-t es b-t szamolja
mikro=mean(x);
nu=mean(y);
Qxy=sum((x-mikro).*(y-nu));
Qx=sum((x-mikro).^2);
Qy=sum((y-nu).^2);
b=Qxy/Qx;
a=nu-b*mikro;

% detrend
yy=a+x.*b;
result=y-yy;

