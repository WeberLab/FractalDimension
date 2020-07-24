function [Hurst] = FDOlga(functional,TR)
%matlabpool close
parpool
% Parameters:
%TR = 2.5; % TR = 2600ms was used
% Load Data
%anatomical = load_nii('1399_anat_ns.nii');
%functional = load_nii('363func_ALIGNED.nii');
[N1,N2,N3,~] = size(functional);

% Fractal Dimension per voxel
Hurst = zeros(N1,N2,N3);
Hurst_fGN = zeros(N1,N2,N3);
Hurst_fBM = zeros(N1,N2,N3);
FD = zeros(N1,N2,N3);
%FDstd = FD;
%method = 'fgm';

parfor i = 1:N1
    tic
    for j = 1:N2
        for k = 1:N3
%             if sum(abs(squeeze(functional(i,j,k,:))))>0
%                 [FD{i,j,k} FD2{i,j,k}] = fractaldim(squeeze(functional.img(i,j,k,:)),0,3);
%             end
            rawBOLD = double(squeeze(functional(i,j,k,:)));
            [Hurst(i,j,k)] = HurstOlga(rawBOLD,TR);
        end
    end
    fprintf('Loop [%g/%g]... %g s\n',i,N1,toc);
end

% FD(isnan(FD)) = 0;
% FD(FD==0) = 1/2;
Hurst(isnan(Hurst)) = 0;
Hurst(Hurst==0) = 1/2;
