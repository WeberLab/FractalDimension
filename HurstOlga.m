function [Hurst] = HurstOlga(rawBOLD, TR)
% This function will analyze an ROI timecourse and output the appropriate
% Hurst exponent according to the procedure of Eke et al. (Eur J Physiol
% (2000) 439:403?415)
% Olga Dona (2016/02/20), Saurabh Shaw(2015), Mohamed A. Warsi(2012), Alex
%Weber(2012), Alya Elzybaak(2011), Evan McNabb(2013), Michael Noseworthy.
% Initialziation:
disp('test')
%Hurst = NaN;

%Hurst_SSC_fGN_Dispersion = NaN;
%Hurst_SSC_fBM_SWV = NaN;
%H = NaN;

abort = false;

% Common parameters:
fs = 1/TR; % Sampling frequency
disp(['sampling frequency = ', num2str(fs)])
n = length(rawBOLD); % Number of timepoints
disp(['number of timepoints = ', num2str(n)])

% Normalizing the time series:
m1 = mean(rawBOLD);
disp(['mean rawBOLD = ', num2str(m1)])
rawBOLD_sub = rawBOLD - m1;

% Multiply each new value by parabolic window
N = length(rawBOLD_sub);
W = zeros(N, 1);
for j = 1:N
    W(j) = 1 - (2*j/(N+1)-1).^2; % parabolic window
end
rawBOLD_sub_double = double(rawBOLD_sub); %convert form int16 to double
signal_pw = rawBOLD_sub_double.*W;

% Matching the ends:
y11 = signal_pw(1); y21 = signal_pw(end);
slope1 = (y21-y11)/(N-1);
y_int1 = y21 - slope1*N;
line = 1:N;
E1 = slope1 * line + y_int1;

% Bridge detrend:
signal_em1 = signal_pw - E1';
range = ceil((N+1) / 2);
freq = [fs * (0 : range-1) / N]';

% plot log (power) vs. log(frequency) --> make sure this is linear over a
% 2-decade range otherwise signal can't be analyzed using fractals
fftSignal1 = fft(signal_em1,N);
fftSignal1 = fftSignal1(1:range); % 1st half of fft since it's symmetric
PSD1 = (abs(fftSignal1).^2)/N;

if rem(N,2)
    PSD1(2:end)=PSD1(2:end)*2;
else
    PSD1(2:end-1)=PSD1(2:end-1)*2;
end
[~, min_index1] = min(abs(freq - 0.08));
[~, min_index2] = min(abs(freq - 0.16));
PSDpart = PSD1(min_index1 : min_index2); % fit only low frequencies for the fit (0.08
%Hz - 0.16 Hz). This is to exculde random white noise at higher frequencies and magnet
%noise below 0.02 Hz

freqpart = freq(min_index1 : min_index2);
logPSD = log10(PSDpart);
logfreq = log10(freqpart);

% Removing all the Inf values:
nu = ~isinf(logPSD) & ~isinf(logfreq);
logPSD_c = logPSD(nu);
logfreq_c = logfreq(nu);
if (~isempty(logPSD_c) && ~isempty(logfreq_c))
    disp('signal is not noise')
    fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'off');
    [fits_result, fit_goodness] = fit(logfreq_c, logPSD_c, 'poly1', fit_opts);
    
    % Beta is the negative of the slope of the fitted line:
    Beta = -1 * fits_result(1);
    disp(['Beta is: ', num2str(Beta)])
    %RSquare_Beta = fit_goodness.rsquare;
else
    disp('signal is noise')
    % Voxel lies outside the brain, disregard it.
    %abort = true;
    %Beta = 0;
    Hurst = NaN;
    %Hurst_SSC_fGN_Dispersion = NaN;
    %Hurst_SSC_fBM_SWV = NaN;
    return
end

% Analysis if signal is fractional Gaussian Noise (fGn)
if((Beta > -1 && Beta < 0.38) && ~abort)
    disp('signal is fGn')
    %Hurst_PSD_fGn = (Beta + 1) / 2; %This method of calculating H (from slope of line)
    %is not as accurate (according to Eke) as doing the dispersional analysis, so use H from
    %dispersional analysis
    % Dispersional analysis to get H for fGn signals:
    maxBins = nextpow2(length(rawBOLD)) - 1;
    signal_2 = rawBOLD(1 : 2^maxBins);
    
    DISP = zeros(maxBins, 1);
    tau = zeros(maxBins, 1);
    
    for i = 1 : maxBins
        m = 2^i;
        signal_binned = reshape(signal_2, [m, (length(signal_2)/m)]);
        
        mean_binned = mean(signal_binned);
        DISP(i) = std(mean_binned);
        tau(i) = m;
        if(DISP(i) == 0)
            DISP(i) = 0.0001;
        end
    end
    
    logDISP = log10(DISP(3:7));
    logtau = log10(tau(3:7));
    
    fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'off');
    [fits_result, fit_goodness] = fit(logtau, logDISP, 'poly1', fit_opts);
    
    Hurst_fGn_Dispersion = fits_result.p1 + 1;
    %RSquare = fit_goodness.rsquare;
    %Hurst = Hurst_PSD_fGn;
    Hurst = Hurst_fGn_Dispersion;
    
    % Analysis if signal is fractional Brownian motion (fBm)
elseif((Beta > 1.04 && Beta < 3) && ~abort)
    disp('signal is fBm')
    %Hurst_PSD_fBm = (Beta - 1) / 2;
    
    %Bridge detrended SWV to get H for fBm signals
    %Recall: signal_em1 is the original bridge detrended data
    maxBins = nextpow2(length(rawBOLD)) - 1;
    signal_2 = signal_em1(1 : 2^maxBins);
    
    SWV = zeros(maxBins, 1);
    tau = zeros(maxBins, 1);
    
    for i = 1 : maxBins
        m = 2^i;
        signal_binned = reshape(signal_2, [m, (length(signal_2)/m)]);
        std_binned = std(signal_binned);
        SWV(i) = mean(std_binned);
        tau(i) = m;
    end
    
    logSWV = log10(SWV(3:7));
    logtau = log10(tau(3:7));
    
    fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'off');
    [fits_result, fit_goodness] = fit(logtau, logSWV, 'poly1', fit_opts);
    
    Hurst_fBm_SWV = fits_result.p1;
    %RSquare = fit_goodness.rsquare;
    %Hurst = Hurst_PSD_fBm;
    Hurst = Hurst_fBm_SWV;
    
    % Signal summation conversion for signals that fall in the non-classifiable region
elseif((Beta >= 0.38 && Beta <= 1.04) && ~abort)
    disp('signal is non-classifiable')
    Y = zeros(size(rawBOLD));
    
    for j = 1:length(Y)
        temp = 0;
        
        for i = 1:j
            temp = temp + rawBOLD(i);
        end
        
        Y(j) = temp;
    end
    
    
    % Run SWV to get Hurst exponent:
    maxBins = nextpow2(length(rawBOLD)) - 1;
    signal_2 = signal_em1(1 : 2^maxBins);
    
    SWV = zeros(maxBins, 1);
    tau = zeros(maxBins, 1);
    
    for i = 1 : maxBins
        m = 2^i;
        signal_binned = reshape(signal_2, [m, (length(signal_2)/m)]);
        
        std_binned = std(signal_binned);
        SWV(i) = mean(std_binned);
        tau(i) = m;
    end
    logSWV = log10(SWV(3:7));
    logtau = log10(tau(3:7));
    
    fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'off');
    [fits_result, fit_goodness] = fit(logtau, logSWV, 'poly1', fit_opts);
    
    Hurst= fits_result.p1;
    %Stats_RSquare = fit_goodness.rsquare;
    disp(['Hurst after summation is: ', num2str(Hurst)])
    
    if(Hurst < 0.8)
        % The signal is an fGn signal, so can do dispersion analysis on it
        % to get final Hurst
        disp('signal is now classified as fGn')
        maxBins = nextpow2(length(rawBOLD)) - 1;
        signal_2 = rawBOLD(1 : 2^maxBins);
        
        DISP = zeros(maxBins, 1);
        tau = zeros(maxBins, 1);
        
        for i = 1 : maxBins
            m = 2^i;
            signal_binned = reshape(signal_2, [m, (length(signal_2)/m)]);
            
            mean_binned = mean(signal_binned);
            DISP(i) = std(mean_binned);
            tau(i) = m;
            if(DISP(i) == 0)
                DISP(i) = 0.0001;
            end
        end
        
        logDISP = log10(DISP(3:7));
        logtau = log10(tau(3:7));
        
        fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'off');
        [fits_result, fit_goodness] = fit(logtau, logDISP, 'poly1', fit_opts);
        
        Hurst_SSC_fGN_Dispersion = fits_result.p1 + 1;
        %RSquare = fit_goodness.rsquare;
        Hurst = Hurst_SSC_fGN_Dispersion;
    elseif(Hurst > 1)
        %The signal is an fBm signal, so can do SWV analysis on it to get
        %final Hurst
        disp('signal is now classified as fBm')
        maxBins = nextpow2(length(rawBOLD)) - 1;
        signal_3 = signal_em1(1 : 2^maxBins);
        
        SWV1 = zeros(maxBins, 1);
        tau1 = zeros(maxBins, 1);
        
        for v = 1 : maxBins
            mm = 2^v;
            signal_bin = reshape(signal_3, [mm, (length(signal_3)/mm)]);
            
            std_bin = std(signal_bin);
            SWV1(i) = mean(std_bin);
            tau1(i) = mm;
        end
        
        logSWV1 = log10(SWV1(3:7));
        logtau1 = log10(tau1(3:7));
        
        fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'off');
        [fits_result, fit_goodness] = fit(logtau1, logSWV1, 'poly1', fit_opts);
        
        Hurst_SSC_fBM_SWV= fits_result.p1;
        %RSquare = fit_goodness.rsquare;
        Hurst = Hurst_SSC_fBM_SWV;
    else
        % If H is not less than 0.8 or larger than 1, then it can't be
        % classified
        Hurst = NaN;
        disp('H is between 0.8 and 1 and cannot be classified')
    end
else
    disp('signal is noise?')
    Hurst = NaN;

end