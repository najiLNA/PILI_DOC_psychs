%This script finds the optimal a values for a model fit at the indvidual subject level to defined
%structural and functional connectivities and an optimal G value .

clear
close all

nReps = 10;

fs = 1/2; %Sampling frequency 1/TR

fc1 = 0.04; % Lower bound for the band-pass filter
fc2 = 0.07; % Upper bound for the band-pass filter


% Group to process
selected_group = 'UWS';

dirDTI = [pwd '\structural_data'];
dirG = [pwd '\Optimal G per subject ' selected_group '\'];


% Order of the subjects
order_MCS = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
   ,'S16','S17','S18','S19','S20','S21','S22','S23', 'S24', 'S25','S26'};

order_UWS = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
   ,'S16','S17','S18','S19','S20'};
order_CNT = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14'...
    ,'S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'...
    'S31','S32','S33','S34','S35'};

%load in indexes of alphas of resting state networks
load([ pwd '\data_ICA_9.mat']);

dirDTI_MCS = [dirDTI '\MCS_structure_mat_files\'];
dirDTI_UWS = [dirDTI '\UWS_structure_mat_files\'];
dirDTI_CNT = [dirDTI '\CNT_structure_mat_files\'];

% Load bold
dataDir = [pwd '\BOLD_data\Data_DoC\' selected_group '\'];
recordings=dir([dataDir '*.mat']);
genericMatrix = struct2cell(load([dataDir '/' recordings(1).name]));
time_series = genericMatrix{1};

dirSave = [pwd '\Optimal a per subject' selected_group];

switch selected_group
    case 'CNT'
        dirDTI = dirDTI_CNT;
        order = order_CNT;
    case 'UWS'
        dirDTI = dirDTI_UWS;
        order = order_UWS;
    case 'MCS'
        dirDTI = dirDTI_MCS;
        order = order_MCS;
end


for nreg=1:length(order)

    % Load G
    load([dirG '/' order{nreg} '.mat']);
    G = optimalG;    

    % First load the subject
    W1.time_series = time_series(nreg,:,:);
  
    W1.name = order{nreg};
    
    % Reorder time series according to symmetrical to match SC
    W1_EO = Deco2AAL(squeeze(W1.time_series)); % Fixed the code to work with our AAL order
    
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
    
    %% HOPF MODEL
    % Load the SC matrix
    % Obtained from the actual DTI
    fprintf('Subject: %s %d\n', recordings(nreg).name(1:7), nreg);
    SC_selected = load([dirDTI '/' order{nreg} '.mat']);
   
    name = [selected_group '_' order{nreg} '_structure'];
    SC = SC_selected.structure;
    SC = SC/max(max(SC))*0.2; % Scale SC matrix to 0.2 
    
    % Parameters for the Hopf Model
    N = length(SC);
    
    TR = 1/fs; % Repetition time (1/fs)
    Tmax = 300; % Time to be simulated in volumes. 
    omega = repmat(2*pi*f_peaks,1,2);  % Changes the frequencies from Hz to radians and copies the array to a second dimension
    omega(:,1) = -omega(:,1); % Change the the sign of the first column to the opposite 
    dt=0.1.*(TR/2); % Compute the time steps from the repetition time to generate the simulated data
    
    sigma = 0.04; % Fixed value of sigma. 
    dsig = sigma*sqrt(dt); 
    Nnew = size(SC,1); % Number of regions (ROIs) in the SC matrix. 90 in the AAL atlas
    isd = find(tril(ones(Nnew),-1)); % This tells us the indexes of the matrix that are in the lower diagonal 
    TTT = Tmax; % Redundant
    
    C = G*SC; % Multiplies the SC matrix by the G value
    
    index = 1;
    lineLength = 0;
    
    grouping_prior = data_ICA_9;
    
    %% Genetic algorithm to obtain optimal a values according to resting state networks
    % We need to define the Hopf model with the input parameter being the 6
    % anatomical priors, which need to be optmized. 
    
    lb = -0.2*ones(1,size(grouping_prior,2)); % Lower bounds for a
    ub = 0.2*ones(1,size(grouping_prior,2)); % Upper bounds for a
    dimx = size(grouping_prior,2); % Number of parameters to optimize
    opts = optimoptions('ga','PlotFcn',{@gaplotbestf,@gaplotstopping,@gaplotbestindiv,@gaplotgenealogy},'Display','off','UseParallel', true, 'UseVectorized', true);
    sumC = repmat(sum(C,2),1,2); % Sums the weighted SC values over on dimension (average connectivity of each ROI, and then creates a new column with the same values. This is needed for the vectorization of the Hopf model
    
    fx_opt = @(a_vector)HopfModel_SSIM_vectorized_eig(C, sumC, a_vector, omega, dsig, Nnew, fc1, fc2, TTT, TR, grouping_prior, FC);
    
    for rep = 1:nReps
        fprintf('Current repetition: %i\n', rep)
        [solution,fval,exitflag,output,population,scores] = ga(fx_opt,dimx,[],[],[],[],lb,ub,[],opts);
        a_vector(rep,:) = solution;
    end
    
    % Save the optimal a values
    mkdir(dirSave)
    save([dirSave recordings_(nreg).name(1:7) '.mat'], 'a_vector')
    
    clearvars W1_EO FC f_peaks_ W1_EO_corr
    
end