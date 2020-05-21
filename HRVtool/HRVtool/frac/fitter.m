function [cs,mcs]=fitter(bs,fs)

i=length(bs);
j=length(fs);
csmax=0;
locmax=0;
while i>2
   i=i-1;
   p = polyfit( log10(bs(1:i)), log10(fs(1:i)), 1 );
   fast(j-i) = p(1);
   uj=polyval(p,log10(bs(1:i)));
   ccs=corrcoef(uj,log10(fs(1:i)));
   cs(j-i,1)=ccs(1,2);
   p = polyfit( log10(bs(i:j)), log10(fs(i:j)), 1 );
   slow(j-i) = p(1);
   uj=polyval(p,log10(bs(i:j)));
   ccs=corrcoef(uj,log10(fs(i:j)));
   cs(j-i,2)=ccs(1,2);
   cs(j-i,3)=cs(j-i,1)+cs(j-i,2);
   if cs(j-i,3)>csmax
      csmax=cs(j-i,3);
      locmax=j-i;
   end
   mcs(1)=locmax;
   mcs(2)=cs(locmax,1);
   mcs(3)=cs(locmax,2);
   mcs(4)=fast(locmax);
   mcs(5)=slow(locmax);
   mcs(6)=bs(locmax);
end



