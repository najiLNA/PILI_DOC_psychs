%This script finds the optimal G values for a model fit to defined
%structural and functional connectivities. 

clear
close all

fs = 1/2.46; % Sampling frequency (1/TR)
fc1 = 0.04; % Lower bound for the band-pass filter
fc2 = 0.07; % Upper bound for the band-pass filter
G_range = 0:0.01:3;
dirSave = [pwd '\Optimal G Propofol\'];
dirSaveFigures = [pwd '\Optimal G Figures Propofol\'];
dirDTI = [pwd '\structural_data'];

% Load bold
dataDir = [pwd '\BOLD_data\Propofol\S2'];
recordings=dir([dataDir '\*.mat']);

for nreg=1:length(recordings)

    genericMatrix = struct2cell(load([dataDir '/' recordings(nreg).name])); 
    W1.time_series = genericMatrix{1};
    W1.name = recordings(nreg).name(1:end-4); 
    
    % Reorder time series according to symmetrical to match SC
    for nGroup=1:size(W1.time_series,1)
        W1_EO(nGroup,:,:) = Deco2AAL(squeeze(W1.time_series(nGroup,:,:))); 
        
        % Demean and detrend
        W1_EO(nGroup,:,:) = sp04_demean(detrend(squeeze(W1_EO(nGroup,:,:)).')).';
    end
    
    % Filter the fMRI time series between 0.04 and 0.07
    filtered_time_series_complete = NaN(size(W1_EO)); % Create and empty array to put the filtered time series
    
    for nGroup=1:size(W1_EO,1)
        
        [filtered_time_series_complete(nGroup,:,:), ~] = sp03_BandpassSimulatedTS(squeeze(W1_EO(nGroup,:,:)), fc1, fc2, 1/fs);
        
        % Get functional connectivity 
        W1_EO_corr(nGroup,:,:)=corr(squeeze(filtered_time_series_complete(nGroup,:,:))');
        
        % Convert to Fisher's z-values
        W1_EO_corr(nGroup,:,:) = atanh(W1_EO_corr(nGroup,:,:));
        
    end
    
    % Get group FC matrix by averaging the FC matrix over all subjects and then
    %  back-transforming it
    W1_EO_corr_mean=tanh(squeeze(nanmean(W1_EO_corr,1)));
    FC = W1_EO_corr_mean;
    FC(isinf(FC))=0; 
    
    
    for nGroup=1:size(W1_EO,1)
        % Obtain peak frequencies from the time series (needed for the Hopf
        % model)
        [TS_emp, Power_Areas, freq, Ku_emp, Phases] = sp01_GetSpectrum(squeeze(W1_EO(nGroup,:,:)), fc1, fc2, 1/fs);
        % Spectral measures
        [maxpow, maxind] = max(Power_Areas);
        f_peaks_(nGroup,:) = freq(maxind)';
    end
    
    f_peaks_ = mean(f_peaks_,1).'; % Average the frequency peaks across subjects
    
    f_peaks = f_peaks_;
    

    % Load the SC matrix
    SC=load([dirDTI '\avg\CNT_avg_structure.mat']);
    fns = fieldnames(SC);
    SC = SC.(fns{1});
    SC = SC/max(max(SC))*0.2; % Scale SC matrix to 0.2 

    %% HOPF MODEL
    
    % Parameters for the Hopf Model
    bp = 0;
    N = length(SC);
    
    a=ones(N,2);
    a(:,1)=bp;a(:,2)=bp; % a (equal to 0 for all regions for now)
    
    TR = 1/fs; % Repetition time (1/fs)
    Tmax = 300; % Time to be simulated in volumes. 
    omega = repmat(2*pi*f_peaks,1,2);  % Changes the frequencies from Hz to radians and copies the array to a second dimension
    omega(:,1) = -omega(:,1); % Change the the sign of the first column to the opposite 
    dt=0.1.*(TR/2); % Compute the time steps from the repetition time to generate the simulated data
    sigma = 0.04;
    dsig = sigma*sqrt(dt);
    Nnew = size(SC,1);
    isd = find(tril(ones(Nnew),-1));
    TTT = Tmax; 
    
    index = 1;
    lineLength = 0;
    
    %% Hopf bifurcation
    for G=G_range % Global coupling parameter (from 0 to 3 in steps of 0.01)
        
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('Current G: %.2f\n',G); % See which G value is being evaluated
        
        parfor sNumber = 1:500 % Simulate the model 500 times
            C = G*SC; % Multiplies the SC matrix by the current WE (G) value (ranging from 0 to 3 in steps of 0.01)
            sumC = repmat(sum(C,2),1,2); % Sums the weighted SC values over one dimension (average connectivity of each ROI, and then creates a new column with the same values. This is needed for the vectorization of the Hopf model
            [FC_sim, ~, ku] = sp00_HopfModel_v3(C, sumC, a, omega, dsig, Nnew, fc1, fc2, TTT, TR); % Hopf model
            FC_sim(isinf(FC_sim))=0;
            % Compute KS distance
            FC_specific_group = squeeze(FC(:,:));
            [h,p,KSdist(index,sNumber)] = kstest2(FC_specific_group(isd),FC_sim(isd));
        end
        
        index = index + 1;
    end
    
    %% Display G and get optimal value
    
    options.x_axis = G_range;
    options.handle     = figure(nreg);
    options.color_area = [128 193 219]./255;    
    options.color_line = [ 52 148 186]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'std';
    
    optimalG = displayKSDistance(KSdist,G_range,options);
    
    title(recordings(nreg).name(1:end-4))
    
    mkdir(dirSaveFigures)
    
    savefig([dirSaveFigures recordings(nreg).name(1:end-4) '.fig'])
    
    % Save the optimal G
    mkdir(dirSave)
    save([dirSave recordings(nreg).name], 'optimalG', 'KSdist')
    
    clearvars W1_EO FC f_peaks_ W1_EO_corr
    
end