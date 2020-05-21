function sd=swvhelp(ts,n,i);   
counter=1;
sd=[];
for j=1:2^i:n-i,
	tsmod=ts(j:j+2^i-1);
%  tsid=bridge(tsmod);
   
   nn=2^i;
   points=(nn-1:-1:0)';
	line=((tsmod(1)-tsmod(nn))*points/(nn-1))+tsmod(nn);
	tsmod=tsmod-line;

%  sd(counter)=std(tsid);
   
   mikro=mean(tsmod);
	sd(counter)=(sum((tsmod-mikro).^2)/(nn-1))^0.5;

   
	tsid=[];
	counter=counter+1;
end;	
