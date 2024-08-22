%This script finds the optimal a values for a model fit at the group level to defined
%structural and functional connectivities and an optimal G value .

clear
close all

nReps = 10;
runNumber = 1;

fs = 1/2.46; % Sampling frequency (1/TR)
fc1 = 0.04; % Lower bound for the band-pass filter
fc2 = 0.07; % Upper bound for the band-pass filter

% Save folder
dirSave = [pwd '\Optimal a Propofol\'];
dirDTI = [pwd '\structural_data\'];
dirG = [pwd '\Optimal G Propofol\'];

% Load bold
dataDir = [pwd '\BOLD_data\Propofol\S2\'];
recordings=dir([dataDir '*.mat']);
% load index file for each resting state network
load([ pwd '\data_ICA_9.mat']);


for nreg=1:length(recordings)

    
    % Load optimal G
    load([dirG '/' recordings(nreg).name]);
    G = optimalG;
    
    genericMatrix = struct2cell(load([dataDir '/' recordings(nreg).name]));
    W1.time_series = genericMatrix{1};
    W1.name = recordings(nreg).name(1:end-4);
    
    % Reorder time series to symmetrical to match SC
    for nGroup=1:size(W1.time_series,1)
        W1_EO(nGroup,:,:) = Deco2AAL(squeeze(W1.time_series(nGroup,:,:))); 
        
        % Demean and detrend
        W1_EO(nGroup,:,:) = sp04_demean(detrend(squeeze(W1_EO(nGroup,:,:)).')).';
    end
    
    % Filter the fMRI time series between 0.04 and 0.07
    filtered_time_series_complete = NaN(size(W1_EO)); 
    
    for nGroup=1:size(W1_EO,1)
        
        [filtered_time_series_complete(nGroup,:,:), ~] = sp03_BandpassSimulatedTS(squeeze(W1_EO(nGroup,:,:)), fc1, fc2, 1/fs);
        W1_EO_corr(nGroup,:,:)=corr(squeeze(filtered_time_series_complete(nGroup,:,:))');
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
    N = length(SC);

    TR = 1/fs; % Repetition time (1/fs)
    Tmax = 10000; % Time to be simulated in volumes. 
    omega = repmat(2*pi*f_peaks,1,2);  % Changes the frequencies from Hz to radians and copies the array to a second dimension
    omega(:,1) = -omega(:,1); % Change the the sign of the first column to the opposite
    dt=0.1.*(TR/2); % Compute the time steps from the repetition time to generate the simulated data
    sigma = 0.04; 
    dsig = sigma*sqrt(dt); 
    Nnew = size(SC,1);
    isd = find(tril(ones(Nnew),-1));
    TTT = Tmax;    
    C = G*SC; % Multiplies the SC matrix by the G value
    index = 1;
    lineLength = 0;
    grouping_prior = data_ICA_9;
  
    % Define how the a parameter works with the anatomical priors
    %a_prior = @(a_vector,grouping) a_vector*grouping.';
    
    %% Genetic algorithm to obtain optimal a values according to anatomical priors
    % We need to define the Hopf model with the input parameter being the 
    % resting state networks, which need to be optmized. 
    
    lb = -0.2*ones(1,size(grouping_prior,2)); % Lower bounds for a
    ub = 0.2*ones(1,size(grouping_prior,2)); % Upper bounds for a
        
    dimx = size(grouping_prior,2); % Number of parameters to optimize
    opts = optimoptions('ga','PlotFcn',{@gaplotbestf,@gaplotstopping,@gaplotbestindiv,@gaplotgenealogy},'Display','off','UseParallel', true, 'UseVectorized', true);
    sumC = repmat(sum(C,2),1,2); % Sums the weighted SC values over on dimension (average connectivity of each ROI, and then creates a new column with the same values. This is needed for the vectorization of the Hopf model
    fx_opt = @(a_vector)HopfModel_SSIM_vectorized_eig(C, sumC, a_vector, omega, dsig, Nnew, fc1, fc2, TTT, TR, grouping_prior, FC);
    
   
    for rep = 1:nReps
        disp('##################################')
        disp(rep)
        disp(recordings(nreg).name)
        disp('##################################')
        
        [solution,fval,exitflag,output,population,scores] = ga(fx_opt,dimx,[],[],[],[],lb,ub,[],opts);
        a_vector(rep,:) = solution;
        
    end
    
    % Save the optimal a values
    mkdir(dirSave)
    save([dirSave recordings(nreg).name], 'a_vector')
    
    clearvars W1_EO FC f_peaks_ W1_EO_corr
    
end