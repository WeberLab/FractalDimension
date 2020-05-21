eval(['load ' '1003_ROI;']);
SLnum=1;
ROIname='righthem';
BW_interp=imresize(BW,'scale', 8, 'method', 'nearest'); % % MAW_Feb2010 used for saving ROI data

ROIbold_mean = BW.*bold_mean(:,:,SLnum);
ROIbold_std = BW.*bold_std(:,:,SLnum);
ROIFDmapPS = BW.*FDmapPS(:,:,SLnum);
ROIFDmapRDall = BW.*FDmapRDall(:,:,SLnum);      % % MAW_Feb2010
ROI_all_interp = BW_interp.*FDmapRDall_interp(:,:,SLnum);    % % MAW_Feb2010
ROI_PS_interp = BW_interp.*FDmapPS_interp(:,:,SLnum);    % % MAW_Feb2010

[qhy1,qhx1,MN]=find(ROIbold_mean);
[qhy2,qhx2,stdev]=find(ROIbold_std);
[qhy5,qhx5,FDps]=find(ROIFDmapPS);
[qhy6,qhx6,FDRDall]=find(ROIFDmapRDall);        % % MAW_Feb2010
[qhy9,qhx9,ROI_all_interp]=find(ROI_all_interp);        % % MAW_Feb2010
[qhy10,qhx10,ROI_PS_interp]=find(ROI_PS_interp);        % % MAW_Feb2010

mean_BOLD = mean(MN); median_BOLD = median(MN); std_BOLD = std(MN); jbtest_BOLD = jbtest(MN,0.01);
mean_stdev = mean(stdev); median_stdev = median(stdev); std_stdev = std(stdev); jbtest_stdev = jbtest(stdev,0.01);
mean_FDRDall = mean(FDRDall); median_FDRDall = median(FDRDall); std_FDRDall = std(FDRDall); jbtest_FDRDall = jbtest(FDRDall,0.01);
mean_FDps = mean(FDps); median_FDps = median(FDps); std_FDps = std(FDps); jbtest_FDps = jbtest(FDps,0.01);
RobustFDStats(:,:,SLnum,3)= BW.*RobustFDStats(:,:,SLnum,3);
RobustFDStats(:,:,SLnum,4)= BW.*RobustFDStats(:,:,SLnum,4);
ROI_size=sum(sum(BW));
sig_fit_all = ROI_size-(sum(sum(sum(RobustFDStats(:,:,SLnum,3)>0.05))));
sig_fit_PS = ROI_size-(sum(sum(sum(RobustFDStats(:,:,SLnum,4)>0.05))));


statsdata=[mean_BOLD,median_BOLD,std_BOLD,jbtest_BOLD;
    mean_stdev,median_stdev,std_stdev,jbtest_stdev;
    mean_FDRDall,median_FDRDall,std_FDRDall,jbtest_FDRDall;
    mean_FDps,median_FDps,std_FDps,jbtest_FDps;
    ROI_size,0,sig_fit_all,sig_fit_PS];

ROIname=char(ROIname);
eval(['save roistats.txt' ' statsdata -ASCII;']);
clear mean_BOLD median_BOLD std_BOLD jbtest_BOLD;
clear mean_stdev median_stdev std_stdev jbtest_stdev;
clear mean_FDRDslow median_FDRDslow std_FDRDslow jbtest_FDRDslow;
clear mean_FDRDfast median_FDRDfast std_FDRDfast jbtest_FDRDfast;
clear mean_FDRDall median_FDRDall std_FDRDall jbtest_FDRDall;
clear mean_FDps median_FDps std_FDps jbtest_FDps;
clear sig_fit_fast sig_fit_slow sig_fit_all sig_fit_PS;