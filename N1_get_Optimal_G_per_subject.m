%This script finds the optimal G values for a model fit to indvidually defined
%structural and functional connectivities. 


clear
close all

fs = 1/2; % Sampling frequency (1/TR)

fc1 = 0.04; % Lower bound for the band-pass filter
fc2 = 0.07; % Upper bound for the band-pass filter

% G range
G_range_top_number = 3;
G_range = 0:0.01:G_range_top_number;

% orderkey
order_MCS = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
   ,'S16','S17','S18','S19','S20','S21','S22','S23', 'S24', 'S25','S26'};

order_UWS = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
   ,'S16','S17','S18','S19','S20'};

order_CNT = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
   ,'S16','S17','S18','S19','S20','S21','S22', 'S23', 'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33', 'S34', 'S35'};

% Save folder
dirSave = [pwd '\Optimal G per subject UWS\'];
dirSaveFigures = [pwd '\Optimal G per subject UWS\'];
dirDTI = [pwd '\structural_data'];

dirDTI_MCS = [dirDTI '\MCS_structure_mat_files\'];
dirDTI_UWS = [dirDTI '\UWS_structure_mat_files\'];
dirDTI_CNT = [dirDTI '\CNT_structure_mat_files\'];


% Load bold
dataDir = [pwd '\BOLD_data\Data_DoC\UWS\'];
recordings=dir([dataDir '*.mat']);

genericMatrix = struct2cell(load([dataDir '\' recordings(1).name]));

for nreg=1:length(order_UWS)
    
    % First load the subject
    W1.time_series = genericMatrix{1};
    W1.time_series = W1.time_series(nreg,:,:);
    W1.name = order_MCS(nreg);
    
    % Reorder time series according to symmetrical to match SC
    W1_EO = Deco2AAL(squeeze(W1.time_series)); 
    
    % Demean and detrend
    W1_EO = sp04_demean(detrend(squeeze(W1_EO).')).';
    
    % Filter the fMRI time series between 0.04 and 0.07
    filtered_time_series_complete = NaN(size(W1_EO));
    
    [filtered_time_series_complete, ~] = sp03_BandpassSimulatedTS(squeeze(W1_EO), fc1, fc2, 1/fs);
    
    % Get functional connectivity by computing the correlation between the
    % time series
    W1_EO_corr = corr(squeeze(filtered_time_series_complete)');

    % Get group FC matrix by averaging the FC matrix over all subjects and then
    %  back-transforming it
    FC = W1_EO_corr;
    
    FC(isinf(FC)) = 0; % Since Fisher's transform turns 1s to Inf, change the diagonal to 0 
    
    
    % Obtain peak frequencies from the time series (needed for the Hopf
    % model)
    [TS_emp, Power_Areas, freq, Ku_emp, Phases] = sp01_GetSpectrum(squeeze(W1_EO), fc1, fc2, 1/fs);
    % Spectral measures
    [maxpow, maxind] = max(Power_Areas);
    f_peaks_ = freq(maxind)'; 
    f_peaks = f_peaks_;

    % Load the SC matrix
    
    
     
    %Assign SC based on diagnosis group
    switch recordings(nreg).name(5:7)
    
        case 'MCS'
            SC = load([dirDTI_MCS  order_MCS{nreg} '.mat']);
        case 'UWS'
            SC = load([dirDTI_UWS  order_UWS{nreg} '.mat']);
        case 'CNT'
            SC = load([dirDTI_CNT  order_CNT{nreg} '.mat']);
    
    end
    
    SC = SC.('structure');
    
    %% HOPF MODEL
   
    SC = SC/max(max(SC))*0.2; % Scale SC matrix to 0.2
    % Parameters for the Hopf Model
    bp = 0;
    N = length(SC);
    
    a=ones(N,2);
    a(:,1)=bp;a(:,2)=bp; % a (equal to 0 for all regions for now)
    
    TR = 1/fs; % Repetition time (1/fs)
    Tmax = 300; 
    omega = repmat(2*pi*f_peaks,1,2);  
    omega(:,1) = -omega(:,1); 
    dt=0.1.*(TR/2); % Compute the time steps from the repetition time to generate the simulated data
    
    sigma = 0.04; 
    dsig = sigma*sqrt(dt); 
    Nnew = size(SC,1); 
    isd = find(tril(ones(Nnew),-1)); % This tells us the indexes of the matrix that are in the lower diagonal (which has all the info since the SC matrix is symmetrical
    TTT = Tmax; 
    
    index = 1;
    lineLength = 0;
    
    recordings(nreg).name
    fprintf('subject: %d ', nreg);

    %% Hopf bifurcation
    for G=0:0.01:G_range_top_number % Global coupling parameter (from 0 to 3 in steps of 0.01)
        
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('Current G: %.2f\n',G); % See which G value is being evaluated
        
        parfor sNumber = 1:100 % Simulate the model 50 times like in the paper
            C = G*SC; % Multiplies the SC matrix by the current WE (G) value (ranging from 0 to 31 in steps of 0.01)
            sumC = repmat(sum(C,2),1,2); % Sums the weighted SC values over one dimension (average connectivity of each ROI, and then creates a new column with the same values. This is needed for the vectorization of the Hopf model
            
            [FC_sim, ~, ku] = sp00_HopfModel(C, sumC, a, omega, dsig, Nnew, fc1, fc2, TTT, TR); % Hopf model
            
          
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
    options.color_area = [128 193 219]./255;    % Blue theme
    options.color_line = [ 52 148 186]./255;
    options.alpha      = 0.5;
    options.line_width = 3;
    options.error      = 'std';
    
    optimalG = displayKSDistance(KSdist,G_range,options);
    
    title(order_MCS(nreg))
    
    mkdir(dirSaveFigures)
    
    savefig([dirSaveFigures order_MCS{nreg} '.fig'])
    
    % Save the optimal G
    mkdir(dirSave)
    save([dirSave order_MCS{nreg}], 'optimalG', 'KSdist')
    
    clearvars W1_EO FC f_peaks_ W1_EO_corr
    
end