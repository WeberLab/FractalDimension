function ekgmhl(action)


Hz=200; 

global BaseHL 	
global FANA		
global SUNA		
global ID		
global YEAR		
global MONTH	
global DAY		

BaseHL=findobj('Tag','base'); 				
FNNSHL=findobj(BaseHL, 'Tag','SubjM');	 	
FOHL=findobj(BaseHL, 'Tag','OpenM'); 		
FSHL=findobj(BaseHL, 'Tag','SaveM'); 		
FXHL=findobj(BaseHL, 'Tag','ExitM'); 		
FANAT=findobj(BaseHL, 'Tag','FamNameT');	
SUNAT=findobj(BaseHL, 'Tag','SurNameT');	
IDT=findobj(BaseHL, 'Tag','IDT');			
BIDAT=findobj(BaseHL, 'Tag','BirthDayT');	
FANA=findobj(BaseHL, 'Tag','FamName');		
SUNA=findobj(BaseHL, 'Tag','SurName');		
ID=findobj(BaseHL, 'Tag','ID');				
YEAR=findobj(BaseHL, 'Tag','Year');			
MONTH=findobj(BaseHL, 'Tag','Month');		
DAY=findobj(BaseHL, 'Tag','Day');			
OKBNS=findobj(BaseHL, 'Tag','OKButns');	
CCBNS=findobj(BaseHL, 'Tag','CancelButns');
DEBNS=findobj(BaseHL, 'Tag','DelButns');	
SLT=findobj(BaseHL, 'Tag','SubjListT');	
SL=findobj(BaseHL, 'Tag','SubjList');		
DEBD=findobj(BaseHL, 'Tag','DelButnd');	
MOBD=findobj(BaseHL, 'Tag','ModButnd');	
OPBD=findobj(BaseHL, 'Tag','OpeButnd');	
DETT=findobj(BaseHL, 'Tag','DetT');			
DET=findobj(BaseHL, 'Tag','Det');			
DLT=findobj(BaseHL, 'Tag','DataListT');	
DL=findobj(BaseHL, 'Tag','DataList');		
AX=findobj(BaseHL, 'Tag','Axes1');			
SX=findobj(BaseHL, 'Tag','XSlid');			
SY=findobj(BaseHL, 'Tag','YSlid');			
XWM=findobj(BaseHL, 'Tag','XWinMinus');	
XWP=findobj(BaseHL, 'Tag','XWinPlus');		
XWT=findobj(BaseHL, 'Tag','XWinT');			
YWM=findobj(BaseHL, 'Tag','YWinMinus');	
YWP=findobj(BaseHL, 'Tag','YWinPlus');		
YWT=findobj(BaseHL, 'Tag','YWinT');			
SB=findobj(BaseHL, 'Tag','SaveButt');		
PB=findobj(BaseHL, 'Tag','ProcButt');		
FB=findobj(BaseHL, 'Tag','FindButt');		
MB=findobj(BaseHL, 'Tag','MaxBox');			
IB=findobj(BaseHL, 'Tag','InsBox');			
DB=findobj(BaseHL, 'Tag','DelBox');			
STT=findobj(BaseHL, 'Tag','StatT');			
FRT=findobj(BaseHL, 'Tag','FracT');			
RME=findobj(BaseHL, 'Tag','RMeanRR');		
RMET=findobj(BaseHL, 'Tag','RMeanRRT');	
RMA=findobj(BaseHL, 'Tag','RMaxRR');		
RMAT=findobj(BaseHL, 'Tag','RMaxRRT');		
RMI=findobj(BaseHL, 'Tag','RMinRR');		
RMIT=findobj(BaseHL, 'Tag','RMinRRT');		
RMF=findobj(BaseHL, 'Tag','RMeanF');		
RMFT=findobj(BaseHL, 'Tag','RMeanFT');		
RSD=findobj(BaseHL, 'Tag','RSDNN');			
RSDT=findobj(BaseHL, 'Tag','RSDNNT');		
RSA=findobj(BaseHL, 'Tag','RSDANN');		
RSAT=findobj(BaseHL, 'Tag','RSDANNT');		
RRM=findobj(BaseHL, 'Tag','RRMS');			
RRMT=findobj(BaseHL, 'Tag','RRMST');		
RNN=findobj(BaseHL, 'Tag','RNN50');			
RNNT=findobj(BaseHL, 'Tag','RNN50T');		
RPN=findobj(BaseHL, 'Tag','RpNN50');		
RPNT=findobj(BaseHL, 'Tag','RpNN50T');		
RPS=findobj(BaseHL, 'Tag','RPSD');			
RPST=findobj(BaseHL, 'Tag','RPSDT');		
RA1=findobj(BaseHL, 'Tag','RDFAA1');		
RA1T=findobj(BaseHL, 'Tag','RDFAA1T');		
RA2=findobj(BaseHL, 'Tag','RDFAA2');		
RA2T=findobj(BaseHL, 'Tag','RDFAA2T');		
RAB=findobj(BaseHL, 'Tag','RDFAB');			
RABT=findobj(BaseHL, 'Tag','RDFABT');		
RCT=findobj(BaseHL, 'Tag','ClassT');		
RDI=findobj(BaseHL, 'Tag','RDisp');			
RDIT=findobj(BaseHL, 'Tag','RDispT');		
RSW=findobj(BaseHL, 'Tag','RSWV');			
RSWT=findobj(BaseHL, 'Tag','RSWVT');		
RRHB=findobj(BaseHL, 'Tag','RRHistB');		
RRTB=findobj(BaseHL, 'Tag','RRTachB');		
RRSB=findobj(BaseHL, 'Tag','RRScatB');		
RRPB=findobj(BaseHL, 'Tag','RRPhasB');		
RRS2B=findobj(BaseHL, 'Tag','RRScatAB');	
RRP2B=findobj(BaseHL, 'Tag','RRPhasAB');	

switch(action)
   
case 'start' 
   D=dir;
   direxists=0;
   for i=1:length(D);
      if strcmp(D(i).name,'data') & D(i).isdir==1
         direxists=1;
      end
   end
   if ~direxists
      mkdir data
   end
   flist=(getfield(what,'mat'));
   datexists=0;
   for i=1:length(flist);
      if strcmp(char(flist(i)),'ekgbase.mat');
         datexists=1;
      end
   end
   if datexists==0
      subj.fn='';
      subj.sn='';
      subj.id='';
      subj.y=[];
      subj.m=[];
      subj.d=[];
      subj.fi={};
      fil=0;
      save('ekgbase.mat','subj','fil');
   end
   

   
%=================================================================
   
case 'subjinit'
   cre([FANAT SUNAT IDT BIDAT FANA SUNA ID YEAR MONTH DAY OKBNS...
         CCBNS DEBNS SLT SL DEBD MOBD OPBD DETT DET DLT DL]); 	
   set(OKBNS,'UserData',0);
   dis([FSHL])   
   
%=================================================================
  
case 'subjend' 
   des([FANAT SUNAT IDT BIDAT FANA SUNA ID YEAR MONTH DAY OKBNS...
     CCBNS DEBNS SLT SL DEBD MOBD OPBD DETT DET DLT DL AX SX SY]);
   emp([FANA SUNA ID YEAR MONTH DAY]);
   set(OKBNS,'UserData',0);
   set(SL,'Value',1);
   
%=================================================================
   
case 'subjhl'
   ena([OKBNS]);
   load ('ekgbase');
   for i=1:size(subj,2)
      z{i,1}=subj(i).fn;
      z{i,2}=subj(i).sn;
      z{i,3}=subj(i).id;
      z{i,4}=num2str(i);
   end
   zz=sortrows(z,[1 2 3]);
   set(SL,'UserData',zz);
   clear z
   z{1}='New subject';
   for i=1:size(zz,1)
      z{i+1}=zz{i,1};
   end
   set(SL,'String',z);
   v=get(SL,'Value');
   if v==1
      if get(OKBNS,'UserData')==2  | get(SLT,'Userdata')~=v
         emp([FANA SUNA ID YEAR MONTH DAY]);
      end
      set(OKBNS,'UserData',0,'String','New');
      ena([FANA SUNA ID YEAR MONTH DAY]);
      dis([DEBNS DEBD MOBD OPBD DETT DET DLT DL]);
      emp([DET]);
   else
      ena([DEBD MOBD OPBD DETT DET DLT DL]);
      set(FANA,'String',subj(eval(zz{v-1,4})).fn, ...
         'UserData',(eval(zz{v-1,4})));
      set(SUNA,'String',subj(eval(zz{v-1,4})).sn);
      set(ID,'String',subj(eval(zz{v-1,4})).id);
      set(YEAR,'String',subj(eval(zz{v-1,4})).y);
      set(MONTH,'String',subj(eval(zz{v-1,4})).m);
      set(DAY,'String',subj(eval(zz{v-1,4})).d);
      if get(OKBNS,'UserData')==0 | get(OKBNS,'UserData')==1 | ...
            get(SLT,'Userdata')~=v
         dis([FANA SUNA ID YEAR MONTH DAY]);
         ena([DEBNS]);
         set(OKBNS,'UserData',1,'String','Change');
      end
   end
   set(SLT,'UserData',v);
   set(DL,'Value',1);
   
%=================================================================
   
case 'subjclear' 
   load('ekgbase');
   suno=(get(FANA,'UserData'));
   for i=1:size(subj(suno).fi,1)
      ff=subj(suno).fi{i,1}; 
      fn='ek000000';
      fz=size(ff,2);
      ffn=strcat(fn(1:(8-fz)),ff,'.dat');
   	cd .\data
      delete (ffn);
      cd ..
   end
   sn=size(subj,2);
   if suno==1
      sn=sn-1;
      if sn==0
         clear subj
         subj.fn='';
         subj.sn='';
         subj.id='';
         subj.y='';
         subj.m='';
         subj.d='';
         subj.fi={};
         sn=0;
         fil=0;
      end
      if sn>0
         for i=1:sn
            s2(i)=subj(i+1);
         end
         clear subj
         subj=s2;
      end
   else
      sn=sn-1;
      for i=1:suno-1
         s2(i)=subj(i);
      end
      for i=suno:sn
         s2(i)=subj(i+1);
      end
      clear subj
      subj=s2;
   end
   save('ekgbase','subj','fil');
   emp([FANA SUNA ID YEAR MONTH DAY]);
   set(OKBNS,'UserData',0,'String','New');
   set(SL,'Value',1);
      
case 'subjmod'
   [dataerror,allok]=subjdatchk;
   if allok 
      load('ekgbase');
      if (get(OKBNS,'UserData'))==2
         suno=(get(FANA,'UserData')); 
         match=0;
         for i=1:size(subj,2) 
            if strcmp(subj(i).fn,get(FANA,'String')) & strcmp(subj(i).sn,get(SUNA,'String')) & strcmp(subj(i).id,get(ID,'String')) & strcmp(subj(i).y,get(YEAR,'String')) & strcmp(subj(i).m,get(MONTH,'String')) & strcmp(subj(i).d,get(DAY,'String'))
               match=1;
               ekgalert('This subject has already been recorded in the database.');
            end
         end
         if match==0; 
            subj(suno).fn=get(FANA,'String');
            subj(suno).sn=get(SUNA,'String');
            subj(suno).id=get(ID,'String');
            subj(suno).y=get(YEAR,'String');
            subj(suno).m=get(MONTH,'String');
            subj(suno).d=get(DAY,'String');
            sn=size(subj,2);
            fil=fil;
            save('ekgbase','subj','fil');
            emp([FANA SUNA ID YEAR MONTH DAY]);
            set(OKBNS,'UserData',0,'String','Record');
         end
         break 
      end
      if (get(OKBNS,'UserData'))==1
         ena([FANA SUNA ID YEAR MONTH DAY]);
         set(OKBNS,'UserData',2,'String','Save');
      end
      if (get(OKBNS,'UserData'))==0 
         match=0;
         for i=1:size(subj,2) 
            if strcmp(subj(i).fn,get(FANA,'String')) & strcmp(subj(i).sn,get(SUNA,'String')) & strcmp(subj(i).id,get(ID,'String')) & strcmp(subj(i).y,get(YEAR,'String')) & strcmp(subj(i).m,get(MONTH,'String')) & strcmp(subj(i).d,get(DAY,'String'))
               match=1;
               ekgalert('This subject has already been recorded in the database.');
            end
         end
         if match==0
            suno=max(size(subj))+1;
            subj(suno).fn=get(FANA,'String');
            subj(suno).sn=get(SUNA,'String');
            subj(suno).id=get(ID,'String');
            subj(suno).y=get(YEAR,'String');
            subj(suno).m=get(MONTH,'String');
            subj(suno).d=get(DAY,'String');
            subj(suno).fi={};
            sn=suno;
            fil=fil;
            save('ekgbase','subj','fil');
            emp([FANA SUNA ID YEAR MONTH DAY]);
         end
      end
   end
      
%=================================================================
   
case 'deloff'
   dis([DEBNS OKBNS]);
   
%=================================================================
   
case 'dathl'
   if (get(SL,'Value'))==1
      emp([DL]);
      break
   end
   suno=get(FANA,'UserData');
   load ('ekgbase');
   z=subj(suno).fi;
   zz=sortrows(z,[1]);
   for i=1:size(zz,1)
      ff=zz{i,1};
      fz=size(zz{i,1},2);
      fn='ek000000';
      ffn=strcat(fn(1:(8-fz)),ff,'.dat');
      zz{i,1}=ffn;
   end
   set(DL,'UserData',zz);
   clear z
   z{1}='New timeseries';
   for i=1:size(zz,1)
      z{i+1,1}=zz{i,1};
   end
   set(DL,'String',z);
   v=get(DL,'Value');
   if v>size(z,1)|v<0
      v=1;
      set(DL,'Value',1);
   end
   if v==1 
      ena([DET MOBD]);
      emp([DET]);
      dis([DEBD OPBD]);
      set(MOBD,'UserData',0,'String','Record');
   else
      set(DET,'String',zz{v-1,2},'UserData',v-1); 
      if get(MOBD,'UserData')==1 | get(DLT,'Userdata')~=v 
         dis([DET]);
         ena([DEBD OPBD]);
         set(MOBD,'UserData',1,'String','Change');
      end
   end
   set(DLT,'UserData',v);
   
%=================================================================
   
case 'datclear' 
   load('ekgbase');
   suno=(get(FANA,'UserData'));
   num=(get(DET,'UserData'));
   zz=subj(suno).fi;
   ff=subj(suno).fi{num,1}; 
   fn='ek000000';
   fz=size(ff,2);
   ffn=strcat(fn(1:(8-fz)),ff,'.dat');
   cd .\data
   delete (ffn);
   cd ..
   zn=size(zz,1);
   if num==1 
      zn=zn-1;
      if zn==0
         subj(suno).fi={};
      end
      if zn>0 
         for i=1:zn
            z{i,1}=zz{i+1,1};
            z{i,2}=zz{i+1,2};
         end
         clear zz
         zz=z;
         subj(suno).fi=zz;
      end
   else 
      zn=zn-1;
      for i=1:num-1
         z{i,1}=zz{i,1};
         z{i,2}=zz{i,2};
         end
      for i=num:zn
         z{i,1}=zz{i+1,1};
         z{i,2}=zz{i+1,2};
         end
      clear zz
      zz=z;
      subj(suno).fi=zz;
   end
   save('ekgbase','subj','fil');
   emp([DET]);
   set(MOBD,'UserData',0,'String','Record');
   set(DL,'Value',1);
      
%=================================================================
   
case 'datmod'
   load('ekgbase');
   if (get(MOBD,'UserData'))==2 
      suno=(get(FANA,'UserData'));
      num=(get(DET,'UserData'));
      subj(suno).fi{num,2}=get(DET,'String');
      sn=size(subj,2);
      fil=fil;
      save('ekgbase','subj','fil');
      emp([DET]);
      set(MOBD,'UserData',0,'String','Record');
      break 
   end
   if (get(MOBD,'UserData'))==1 
      ena([DET]);
      set(MOBD,'UserData',2,'String','Save');
   end
   if (get(MOBD,'UserData'))==0 
      suno=(get(FANA,'UserData'));
      fil=fil+1;
      sz=size(subj(suno).fi,1)+1;
      subj(suno).fi{sz,1}=num2str(fil);
      subj(suno).fi{sz,2}=get(DET,'String');
      sn=sn;
      ff=num2str(fil);
      fn='ek000000';
      fz=size(ff,2);
      ffn=strcat(fn(1:(8-fz)),ff,'.dat');
      [FileName, PathName] = uigetfile('*.mat', 'Select File');
      if FileName~=0 
         load (strcat(PathName,FileName));
         if ~isempty(TS) 
            sts=size(TS);
            if sts(1)<sts(2)
               TS=TS';
            end
            TSD=zeros(max(size(TS)),2);
            TSD(:,1)=TS;
            cd .\data
            save (ffn,'TSD');
            cd ..
            clear TS TSD
            save('ekgbase','subj','fil');
         end
         
      end
   end
   
%=================================================================
   
case 'workinit' 
   load('ekgbase');
   suno=(get(FANA,'UserData'));
   num=(get(DET,'UserData'));
   ff=subj(suno).fi{num,1};
   fn='ek000000';
   fz=size(ff,2);
   ffn=strcat(fn(1:(8-fz)),ff,'.dat');
   set(BaseHL,'Name',['ECG analysis - ' ffn]);
   cd .\data
   load (ffn,'-mat');
   cd ..
   length=size(TSD);
   set(SX,'UserData',TSD);
   des([FANAT SUNAT IDT BIDAT FANA SUNA ID YEAR MONTH DAY OKBNS...
      CCBNS DEBNS SLT SL DEBD MOBD OPBD DETT DET DLT DL]);
   emp([FANA SUNA ID YEAR MONTH DAY]);
   set(OKBNS,'UserData',0);
   set(SL,'Value',1);
   cre([AX SX SY XWM XWP YWM YWP XWT YWT SB PB FB MB IB DB]);
   dis([MB]);
   set(XWP,'UserData',3); 
   set(YWP,'UserData',20);
   MaxX=size(TSD,1)-get(XWP,'UserData')*Hz;
   set(SX,'Min',1);
   set(SX,'Max',MaxX);
   set(SX,'Value',1);
   set(SX,'SliderStep',[(1/MaxX)*3*Hz*0.75 (1/MaxX)*3*Hz*2]);
   set(SY,'Min',-100);
   set(SY,'Max',100);
   set(SY,'Value',0);
   set(SY,'SliderStep',[0.01 0.05]);
   
%=================================================================
   
case 'xdec'
   s=get(XWP,'UserData'); 
   if s>100
      s=s-100;
   else
      if s>50
         s=s-10;
      else
         if s>10
            s=s-5;
         else
            if s>1
               s=s-1;
            else
               s=1;
            end
         end
      end
   end
   set(XWP,'UserData',s);
   MaxX=size(get(SX,'UserData'),1)-s*Hz;
   set(SX,'SliderStep',[(1/MaxX)*s*Hz*0.75 (1/MaxX)*s*Hz*2]);
   clear s
   
%=================================================================
   
case 'xinc'
   s=get(XWP,'UserData');
   if s<10				
      s=s+1;
   else
      if s<50
         s=s+5;
      else
         if s<100
            s=s+10;
         else
            if s<3000
               s=s+100;
            else
               s=3000;
            end
         end
      end
   end
   set(XWP,'UserData',s);
   MaxX=size(get(SX,'UserData'),1)-s*Hz;
   set(SX,'SliderStep',[(1/MaxX)*s*Hz*0.75 (1/MaxX)*s*Hz*2]);
   clear s
   
%=================================================================
   
case 'ydec'
   s=get(YWP,'UserData');
   if s>100
      s=s-100;
   else
      if s>50
         s=s-10;
      else
         if s>10
            s=s-5;
         else
            if s>1
               s=s-1;
            else
               s=1;
            end
         end
      end
   end
   set(YWP,'UserData',s);
   set(SY,'Min',-(2*s));
   set(SY,'Max',2*s);
   clear s
   
%=================================================================
   
case 'yinc'
   s=get(YWP,'UserData');
   if s<10				
      s=s+1;
   else
      if s<50
         s=s+5;
      else
         if s<100
            s=s+10;
         else
            if s<3000
               s=s+100;
            else
               s=3000;
            end
         end
      end
   end
   set(YWP,'UserData',s);
   set(SY,'Min',-(2*s));
   set(SY,'Max',2*s);
   clear s
   
%=================================================================
   
case 'workslide'
   TSD=get(SX,'UserData');
   l=size(TSD,1);
   sizet=size(TSD);
   xs=get(XWP,'UserData');
   ys=get(YWP,'UserData');
   MaxX=l-xs*Hz;	
   Xp=floor(get(SX,'Value'));
   Yp=floor(get(SY,'Value'));
   MinY=get(SY,'Min');
   MaxY=get(SY,'Max');
   if Xp>MaxX
      Xp=MaxX;
      set(SX,'Value',Xp);
   end
   if Yp>MaxY
      Yp=MaxY;
      set(SY,'Value',Yp);
   end
   if Yp<MinY
      Yp=MinY;
      set(SY,'Value',Yp);
   end
   set(SX,'Max',MaxX);
   xx=(Xp:(Xp+xs*Hz));
   yy=TSD(xx,1);
   z=TSD(xx,2);
   xx=xx/Hz;
   plot(xx,yy);
   zz=find(z==1);
   hold on
   if ~isempty(zz) 
      for i=1:size(zz,1)
         plot([xx(zz(i)) xx(zz(i))], [-ys+Yp ys+Yp],'g-');
      end
   end
   hold off
   zz=find(z==2);
   hold on
   if ~isempty(zz) 
      for i=1:size(zz,1)
         plot([xx(zz(i)) xx(zz(i))], [-ys+Yp ys+Yp],'r-');
      end
   end
   hold off
   axis([Xp/Hz (Xp+xs*Hz)/Hz -ys+Yp ys+Yp]);
   set(gca,'XAxisLocation','top');
   set(gca,'ButtonDownFcn','ekgmhl buttdown, ekgmhl workslide');
      
%=================================================================

case 'buttdown'
   coord=get(gca,'Currentpoint');
   x=floor(coord(1,1)*Hz);
   TSD=get(SX,'UserData');
   del=get(DB,'Value');
   if del>0
      TSD(floor(max(1,(x-0.1*Hz))):floor(min((x+0.1*Hz),max(size(TSD)))),2)=0;
   end
   if get(IB,'Value')>0
      TSD(x,2)=2;
      if get(MB,'Value')>0
         if (x>1)
            while TSD(x,1)<=TSD(x-1,1) 
               TSD(x,2)=0;
               TSD(x-1,2)=2;
               x=x-1;
               if x==1
                  break
               end
            end
         end
         if (x>1)&(x<max(size(TSD)))
            while TSD(x,1)<TSD(x+1,1) 
               TSD(x,2)=0;
               TSD(x+1,2)=2;
               x=x+1;
               if x==max(size(TSD))
                  break
               end
            end
         end
      end
   end
   set(SX,'UserData',TSD);
   set(gca,'ButtonDownFcn','ekgmhl buttdown, ekgmhl workslide');
   
%=================================================================
   
case 'destall' 
 des([FANAT SUNAT IDT BIDAT FANA SUNA ID YEAR MONTH DAY OKBNS ...
       CCBNS DEBNS SLT SL DEBD MOBD OPBD DETT DET DLT DL])
 des([gca SX SY XWM XWP YWM YWP XWT YWT SB PB FB MB IB DB STT ...
       FRT RME RMET RMA RMAT RMI RMIT RMF RMFT RSD RSDT RSA RSAT])
 des([RRM RRMT RNN RNNT RPN RPNT RPS RPST RA1 RA1T RA2 RA2T RAB ...
    RABT RCT RDI RDIT RSW RSWT RRHB RRTB RRSB RRPB RRS2B RRP2B])
 dis([FSHL])
 set(BaseHL,'Name','ECG analysis');
 cla
   
%=================================================================
   
case 'rfind'  
   TSD=get(SX,'UserData');
   x=TSD(:,1);
   TSD(:,2)=0;
   y=zeros(size(x,1),1);
   y(5:size(x,1))=x(5:size(x,1))-x(1:size(x,1)-4);
   y1=zeros(size(x,1),1);
   y1(9:size(x,1))=y(9:size(x,1))+4*y(8:size(x,1)-1)+6*y(7:size(x,1)-2)+4*y(6:size(x,1)-3)+y(5:size(x,1)-4);
   yc=find(y1>21);
   WB=waitbar(0,'Step 1: QRS-complexes, please wait...');
   set(WB,'Units','pixels');
   pos=get(BaseHL,'Position');
   posw=get(WB,'Position');
	xb=pos(1)+(pos(3)-posw(3))/2;
	yb=pos(2)+(pos(4)-posw(4))/2;
	pos=[xb yb posw(3) posw(4)];
	set(WB,'Position',pos);
   for i=1:max(size(yc))
      waitbar(i/max(size(yc)));
      TSD(yc(i),2)=0;
      if (yc(i)+32)<max(size(y1))
         if ~isempty(find(y1(yc(i):yc(i)+32)<-21))
            yf1=find(y1(yc(i):yc(i)+32)<-21);
            yf1=yf1(1);
            TSD(yc(i),2)=1;
            if ~isempty(find(y1(yc(i)+yf1:yc(i)+32)>21))
               yf2=find(y1(yc(i)+yf1:yc(i)+32)>21);
               yf2=yf2(1);
               TSD(yc(i),2)=1;
               if ~isempty(find(y1(yc(i)+yf1+yf2:yc(i)+32)>21))
                  yf3=find(y1(yc(i)+yf1+yf2:yc(i)+32)>21);
                  yf3=yf3(1);
                  TSD(yc(i),2)=1;
                  if ~isempty(find(y1(yc(i)+yf1+yf2+yf3:yc(i)+32)<-21))
                     TSD(yc(i),2)=0;
                  end
               end
            end
         end
      end
   end
   close(WB);
   yc=find(TSD(:,2)>0);
   for i=1:max(size(yc))
      TSD(yc(i)+1:(yc(i)+25),2)=0;
   end
   yc=find(TSD(:,2)>0);
   WB=waitbar(0,'Step 2: Rmax, please wait...');
   set(WB,'Units','pixels');
   pos=get(BaseHL,'Position');
   posw=get(WB,'Position');
	xb=pos(1)+(pos(3)-posw(3))/2;
	yb=pos(2)+(pos(4)-posw(4))/2;
	pos=[xb yb posw(3) posw(4)];
   set(WB,'Position',pos);
   sizz=max(size(TSD));
   for i=1:max(size(yc))
      waitbar(i/max(size(yc)));
      z=yc(i);
      if (z>1)
         while TSD(z,1)<=TSD(z-1,1) 
            TSD(z,2)=0;
            TSD(z-1,2)=1;
            z=z-1;
            if z==1
               break
            end
         end
      end
      if (z<sizz)
         while TSD(z,1)<TSD(z+1,1) 
            TSD(z,2)=0;
            TSD(z+1,2)=1;
            z=z+1;
            if z==sizz
               break
            end
         end
      end
	end
   close(WB);
   set(SX,'UserData',TSD);
   
%=================================================================
   
case 'delbox'
   if get(DB,'Value')>0 
      set(IB,'Value',0);
      dis(MB);
   end
   
%=================================================================
   
case 'insbox' 
   if get(IB,'Value')>0 
      set(DB,'Value',0);
      set(MB,'Value',1);
      ena(MB);
   else
      dis(MB); 
   end
   
%=================================================================
   
case 'worksave'  
   load('ekgbase');
   suno=(get(FANA,'UserData'));
   num=(get(DET,'UserData'));
   TSD=get(SX,'UserData');
   ff=subj(suno).fi{num,1};
   fn='ek000000';
   fz=size(ff,2);
   ffn=strcat(fn(1:(8-fz)),ff,'.dat');
   cd .\data
   save (ffn,'TSD');
   cd ..
   
%=================================================================
   
case 'process'
   TSD=get(SX,'UserData');
   TSS=TSD(:,2);
   TSF(:,1)=find(TSS>0);
   TSF(:,2)=0;
   TSF2=find(TSS>1);
   if ~isempty(TSF2)
      for i=1:max(size(TSF2))
         TSF(find(TSF(:,1)==TSF2(i)),2)=1;
      end
   end
   TS(:,1)=diff(TSF(:,1));
   TS(:,2)=0;
   TS((find(TSF(:,2)>0)-1),2)=1;
   clear TSD TSS TSF TSF2
   count=0;
   TS2=[];
   i=2;
   WB=waitbar(0,'Seeking appropriate data, please wait...');
   set(WB,'Units','pixels');
   pos=get(BaseHL,'Position');
   posw=get(WB,'Position');
	xb=pos(1)+(pos(3)-posw(3))/2;
	yb=pos(2)+(pos(4)-posw(4))/2;
	pos=[xb yb posw(3) posw(4)];
	set(WB,'Position',pos);
   while i<max(size(TS))
      waitbar(i/max(size(TS)));
      if (TS(i,1)<TS(i-1,1)*(2)) & (TS(i,1)>TS(i-1,1)*(0.5))
         count=count+1;
         TS2(count,1)=TS(i,1);
         if TS(i,2)>0
            mod(count)=1;
         else
            mod(count)=0;
         end
      else
         if count<8192
            count=0;
            mod=[];
            TS2=[];
         else
            i=max(size(TS));
         end
      end
      i=i+1;
   end
   close(WB);
   TS=TS2';
   clear TS2;
   TS=TS/Hz;
   if max(size(TS))<8192
      ekgalert('The available data is too short for analysis');
      set (RRHB,'UserData',[]);
      break
   end
   lTS=max(size(TS));
   TS(1:8192)=TS(1+floor((lTS-8192)/2):8192+floor((lTS-8192)/2));
   mod(1:8192)=mod(1+floor((lTS-8192)/2):8192+floor((lTS-8192)/2));
   clear TSS
   mod=max(size(find(mod>0)));
   if (mod/8192)>0.05
      ekgalert('Number of modifications exceed limit.');
   end
   set (RRHB,'UserData',TS);
   set (SX,'Userdata',0);
   
%=================================================================
   
case 'results'
   TS=get(RRHB,'UserData');
   if isempty(TS)
      break
   end
   if max(size(TS))<8192
      ekgalert('The available data is too short for analysis');
      break
   end
   cla
   des([gca SX SY XWM XWP YWM YWP XWT YWT SB PB FB MB IB DB])
   cre([STT FRT RME RMET RMA RMAT RMI RMIT RMF RMFT RSD ...
         RSDT RSA RSAT RRM RRMT RNN RNNT RPN RPNT RPS RPST RA1])
   cre([RA1T RA2 RA2T RAB RABT RCT RDI RDIT RSW RSWT RRHB ...
         RRTB RRSB RRPB RRS2B RRP2B])
   emp([RME RMA RMI RMF RSD RSA RRM RNN RPN RPS RA1 RA2 RAB ...
         RDI RSW ])
   set(RCT,'String','Classification in progress...');
   ena([FSHL])
   set(RME,'String',[num2str(mean(TS)*1000) ' ms']);
   set(RMI,'String',[num2str(min(TS)*1000) ' ms']);
   set(RMA,'String',[num2str(max(TS)*1000) ' ms']);
   set(RMF,'String',[num2str((max(size(TS))/sum(TS))*60) ' bpm']);
   set(RSD,'String',[num2str(std(TS)*1000) ' ms']);
   win=400;
   i=0;
   WB=waitbar(0,'Determining SDANN, please wait...');
   set(WB,'Units','pixels');
   pos=get(BaseHL,'Position');
   posw=get(WB,'Position');
	xb=pos(1)+(pos(3)-posw(3))/2;
	yb=pos(2)+(pos(4)-posw(4))/2;
	pos=[xb yb posw(3) posw(4)];
	set(WB,'Position',pos);
   while (i+win)<max(size(TS))
      i=i+1;
      waitbar(i/(max(size(TS))-win))
      av(i)=mean(TS(i:(i+win)));
   end
   close (WB)
   set(RSA,'String',[num2str(std(av)*1000) ' ms']);
   clear av win i
   set(RRM,'String',[num2str(sqrt(mean(diff(TS).^2))*1000) ' ms']);
   count=max(size(find(abs(diff(TS))>=0.05)));
   set(RNN,'String',num2str(count));
   set(RPN,'String',strcat(num2str(floor(count/max(size(TS))*10000)/100),' %'));
   l=floor(log2(max(size(TS))));
   l2=2^l;
   TS=TS(1:l2);
   result=spec(bridge(window(TS')),l,200,0);
   beta=result(6);
   set(RPS,'String',num2str(result(6)));
   result=dfa(TS',64,16);
   set(RA1,'String',num2str(result(5)));
   set(RA2,'String',num2str(result(1)));
   set(RAB,'String',num2str(result(2)));
   result=disper(TS',l);
   set(RDI,'String',num2str(result(1)));
   sw=bdswv(TS',l);
   set(RSW,'String',num2str(sw(1)));
   set (RCT,'String','Classification NA');
   dis ([RDI RDIT RSW RSWT]);
   fGn=0;
   fBm=0;
   if beta<0.32
      fGn=1
   end
   if beta>1.16
      fBm=1;
   end
   TSS=cumsum(TS);
   sws=bdswv(TSS',l);
   if ~fBm & ~fGn
      if sws(1)<0.84
         fGn=1;
      end
      if sws(1)>1
         fBm=1;
      end
   end
   if fGn==1
      set (RCT,'String','The signal is fGn');
      ena ([RDI RDIT]);
   end
   if fBm==1
      set (RCT,'String','The signal is fBm');
      ena ([RSW RSWT]);
   end
   
%=================================================================
         
case 'opents'
   [FileName, PathName] = uigetfile('*.mat', 'Select File');
   if FileName~=0 
      load (strcat(PathName,FileName));
      if ~isempty(TS)
         set(RRHB,'UserData',TS);
         set(BaseHL,'Name',['External timeseries analysis- ' FileName]);
         ena(FSHL);
      end
      
   end
      
%=================================================================
      
case 'savets'
   [FileName, PathName] = uiputfile('*.mat', 'Save timeseries as');
   if FileName~=0 
   	TS=get(RRHB,'UserData');
      save(strcat(PathName,FileName),'TS');
   end
   
%=================================================================
   
case 'RRHistB'
   pos=get(BaseHL,'Position');
   TS=get(RRHB,'UserData');
   mint=(floor(min(TS)*100))/100;
   maxt=(floor(max(TS)*100))/100;
   h=figure;
   figure(h);
   set(gcf,'Position',pos,'NumberTitle','off','Name',[get(BaseHL,'Name') ' - Histogram']);
   hist(TS,(floor(min(TS)*100))/100:0.01:(floor(max(TS)*100))/100);
     
%=================================================================
   
case 'RRTachB'
   pos=get(BaseHL,'Position');
   TS=get(RRHB,'UserData');
   h=figure;
   figure(h);
   set(gcf,'Position',pos,'NumberTitle','off','Name',[get(BaseHL,'Name') ' - Tachogram']);
   bar(TS);
   
%=================================================================
   
case 'RRScatB'
   pos=get(BaseHL,'Position');
   TS=get(RRHB,'UserData');
   TSY=TS(2:max(size(TS)));
   TSX=TS(1:max(size(TS))-1);
   h=figure;
   figure(h);
   set(gcf,'Position',pos,'NumberTitle','off','Name',[get(BaseHL,'Name') ' - Scatterplot']);
   hold on;
   axis([(floor(min(TSX)*10)/10) (floor(max(TSX)*10)/10)+0.1 (floor(min(TSY)*10)/10) (floor(max(TSY)*10)/10)+0.1]);
   plot(TSX,TSY,'b:');
   plot(TSX,TSY,'r.');
   hold off
   
%=================================================================
   
case 'RRPhasB'
   pos=get(BaseHL,'Position');
   TS=get(RRHB,'UserData');
   TSY=TS(2:max(size(TS)));
   TSX=TS(1:max(size(TS))-1);
   t=1:max(size(TS)-1);
   ts=1:0.1:max(size(TS))-1;
   xs=spline(t,TSX,ts);
   ys=spline(t,TSY,ts);
   h=figure;
   figure(h);
   set(gcf,'Position',pos,'NumberTitle','off','Name',[get(BaseHL,'Name') ' - Phaseplot']);
   hold on;
   axis([(floor(min(TSX)*10)/10) (floor(max(TSX)*10)/10)+0.1 (floor(min(TSY)*10)/10) (floor(max(TSY)*10)/10)+0.1]);
   plot(xs,ys,'b-');
   plot(TSX,TSY,'r.');
   hold off
   
%=================================================================
   
case 'RRScatAB'
   pos=get(BaseHL,'Position');
   TS=get(RRHB,'UserData');
   TSY=TS(2:max(size(TS)));
   TSX=TS(1:max(size(TS))-1);
   h=figure;
   figure(h);
   set(gcf,'Position',pos,'NumberTitle','off','Name',[get(BaseHL,'Name') ' - Scatterplot animation']);
   hold on;
   axis([(floor(min(TSX)*10)/10) (floor(max(TSX)*10)/10)+0.1 (floor(min(TSY)*10)/10) (floor(max(TSY)*10)/10)+0.1]);
   comet(TSX,TSY);
   plot(TSX,TSY,'b:');
   plot(TSX,TSY,'r.');
   hold off
   
%=================================================================
   
case 'RRPhasAB'
   pos=get(BaseHL,'Position');
   TS=get(RRHB,'UserData');
   TSY=TS(2:max(size(TS)));
   TSX=TS(1:max(size(TS))-1);
   t=1:max(size(TS)-1);
   ts=1:0.1:max(size(TS))-1;
   xs=spline(t,TSX,ts);
   ys=spline(t,TSY,ts);
   h=figure;
   figure(h);
   set(gcf,'Position',pos,'NumberTitle','off','Name',[get(BaseHL,'Name') ' - Phaseplot animation']);
   hold on;
   axis([(floor(min(xs)*10)/10) (floor(max(xs)*10)/10)+0.1 (floor(min(ys)*10)/10) (floor(max(ys)*10)/10)+0.1]);
   comet(xs,ys);
   plot(xs,ys,'b-');
   plot(TSX,TSY,'r.');
   hold off
      
end

%=================================================================

function des(n) 
for i=1:length(n)
   set(n(i),'Visible','off');
end

%=================================================================

function cre(n) 
for i=1:length(n)
   set(n(i),'Visible','on');
end

%=================================================================

function dis(n) 
for i=1:length(n)
   set(n(i),'Enable','off');
end

%=================================================================

function ena(n) 
for i=1:length(n)
   set(n(i),'Enable','on');
end

%=================================================================

function emp(n) 
for i=1:length(n)
   set(n(i),'String','');
end

%=================================================================

function fig = ekgalert(a) 

global BaseHL
hib=msgbox(a,'Error!','error');
set(hib,'Units','pixels');
pos=get(BaseHL,'Position');
posw=get(hib,'Position');
xb=pos(1)+(pos(3)-posw(3))/2;
yb=pos(2)+(pos(4)-posw(4))/2;
pos=[xb yb posw(3) posw(4)];
set(hib,'Position',pos);



%=================================================================

function [dataerror,allok]=subjdatchk
global BaseHL 	
global FANA		
global SUNA		
global ID		
global YEAR		
global MONTH	
global DAY		

allok=1;
dataerror=0;
if length(get(FANA,'String'))==0
   allok=0;
   dataerror=1;
   ekgalert('Family name is missing!');
   break
end
if length(get(SUNA,'String'))==0 
   allok=0;
   dataerror=1;
   ekgalert('First name is missing!');
   break
end
if length(get(ID,'String'))==0
   allok=0;
   dataerror=1;
    ekgalert('ID is missing!');
    break
end
yearok=1;
if length(get(YEAR,'String'))==0 
   allok=0;
   dataerror=1;
   yearok=0;
   ekgalert('Year is missing!');
   break
else
   if length(get(YEAR,'String'))>0 & ((eval(get(YEAR,'String'))<1900)) & dataerror==0
      allok=0;
      dataerror=1;
      yearok=0;
      ekgalert('Year is out of range!');
      break
   end
end
monthok=1;
if length(get(MONTH,'String'))==0
   allok=0;
   dataerror=1;
   monthok=0;
   ekgalert('Month is missing!');
   break
else
   if length(get(MONTH,'String'))>0 & ((eval(get(MONTH,'String'))<1) | (eval(get(MONTH,'String'))>12)) & dataerror==0
      allok=0;
      dataerror=1;
      monthok=0;
      ekgalert('Month is out of range!');
      break
   end
end
if length(get(DAY,'String'))==0 
   allok=0;
   dataerror=1;
   ekgalert('Day is missing!');
   break
else
   if monthok & yearok
      if length(get(DAY,'String'))>0 & ((eval(get(DAY,'String'))<1) | (eval(get(DAY,'String'))>eomday((eval(get(YEAR,'String'))),(eval(get(MONTH,'String'))))))&dataerror==0
         allok=0;
         dataerror=1;
         ekgalert('Day is out of range!');
         break
      end
   end
end
d=datevec(now);
daterr=0;
if eval(get(YEAR,'String'))>d(1)
   allok=0;
   dataerror=1;
   daterr=1;
end
if eval(get(YEAR,'String'))==d(1) & eval(get(MONTH,'String'))>d(2)
   allok=0;
   dataerror=1;
   daterr=1;
end
if eval(get(YEAR,'String'))==d(1) & eval(get(MONTH,'String'))==d(2) & eval(get(DAY,'String'))>d(3)
   allok=0;
   dataerror=1;
   daterr=1;
end
if daterr
   ekgalert('Future date!');
end
