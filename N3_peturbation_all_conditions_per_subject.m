% This script uses the optimal G and a values to build models fit at the individual subject level to
% defined strucutral and functional connectivities for each condition and then calculates PILI
% through perturbing each node of the AAL.

clear
close all

fs = 1/2; % Sampling frequency (1/TR)

fc1 = 0.04; % Lower bound for the band-pass filter
fc2 = 0.07; % Upper bound for the band-pass filter
nSim = 100; % Number of times to simulate the perturbation
nFig = 1;


% Group to process
selected_group = 'UWS';

dirDTI = [pwd '\structural_data'];
dirG = [pwd '\Optimal G per subject ' selected_group '\'];
dirA = [pwd '\Optimal a per subject ' selected_group '\'];

dirDTI_MCS = [dirDTI '\MCS_structure_mat_files\'];
dirDTI_UWS = [dirDTI '\UWS_structure_mat_files\'];
dirDTI_CNT = [dirDTI '\CNT_structure_mat_files\'];

recordings_CNT=dir([dirDTI_CNT '*.mat']);
recordings_MCS=dir([dirDTI_MCS '*.mat']);
recordings_UWS=dir([dirDTI_UWS '*.mat']);

%Drugs to simulate
drugSim = {'Baseline','LSD','Propofol_S2','PSILO_maas', 'PSILO_imperial'};

% Load bold
dataDir = [pwd '\BOLD_data\Data_DoC\' selected_group '\'];
recordings=dir([dataDir '*.mat']);
genericMatrix = struct2cell(load([dataDir '/' recordings(1).name]));
time_series = genericMatrix{1};

grouping = load([ pwd '\data_ICA_9.mat']);
grouping.labels = {'ATTENTION', 'AUD', 'DMN', 'FP', 'HPC', 'PRECUNEUS', 'SM', 'THAL','VIS'};

% Order of the subjects
order_MCS = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
    ,'S16','S17','S18','S19','S20','S21','S22','S23', 'S24', 'S25','S26'};

order_UWS = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
    ,'S16','S17','S18','S19','S20'};
order_CNT = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15' ...
    ,'S16','S17','S18','S19','S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33', 'S34', 'S35'};


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

dirSave = [pwd '\perturbation per subject ' selected_group];
dirSaveFigures = [pwd '\perturbation per subject ' selected_group];


% The delta g and alphas below were calculated manually from the outputs from N2 optima
% alphas scrpt
for drug=1:length(drugSim)
    drugSim{drug}
    switch drugSim{drug}
        case 'Baseline'
            delta_G = 0;
            delta_a = [0 0 0 0 0 0 0 0 0];
        case 'LSD'
            delta_G = 0.04;
            delta_a = [0.0006   -0.0065    0.0430   -0.0016   -0.0024    0.0121   -0.0107    0.0125   -0.0518];
        case 'PSILO_2'
            delta_G = -0.03;
            delta_a =  [-0.0412    0.0927    0.0383   -0.0381    0.0025   -0.0298   -0.0073    0.0218   -0.0936];
        case 'Ket_S1'
            delta_G = -0.22;
            delta_a = [0.0008   -0.0041   -0.0013    0.1207   -0.0188    0.0121    0.0047   -0.0080   -0.0447];
        case 'Ket_S2'
            delta_G = -0.35;
            delta_a = [ -0.0281    0.0596   -0.0041    0.1206   -0.0741    0.0077    0.1639   -0.1797   -0.0112];
        case 'Propofol_S1'
            delta_G = -0.07;
            delta_a = [0.0160    0.0020    0.0436    0.1101    0.0006   -0.0014   -0.0167   -0.0002   -0.1158];
        case 'Propofol_S2'
            delta_G = -0.21;
            delta_a = [ -0.0215   -0.0317   -0.0548   -0.0303   -0.0141    0.0620    0.1166   -0.1866    0.0738];
        case 'PSILO_imperial'
            delta_G = -0.02;
            delta_a = [-0.0459    0.0903    0.0237   -0.0228    0.0132   -0.0074    0.0103    0.0008   -0.0930];
        case 'PSILO_maas'
            delta_G = -0.06;
            delta_a = [ -0.0364    0.0951    0.0528   -0.0534   -0.0083   -0.0522   -0.0250    0.0428   -0.0942];

    end



    for nreg=1:length(order)
   
        % First load the subject
        W1.time_series = (time_series(nreg,:,:));
        W1.name = [recordings(nreg).name(1:end-4) '_' order{nreg}];

        % Load G
        load([dirG order{nreg} '.mat']);
        G = optimalG;
        G_perturb = optimalG + delta_G;

        % Reorder time series to symmetrical to match SC
        for nGroup=1:size(W1.time_series,1)
            W1_EO(nGroup,:,:) = Deco2AAL(squeeze(W1.time_series(nGroup,:,:))); 

            % Demean and detrend    
            W1_EO(nGroup,:,:) = sp04_demean(detrend(squeeze(W1_EO(nGroup,:,:)).')).';
        end

        % Filter the fMRI time series between 0.04 and 0.07
        filtered_time_series_complete = NaN(size(W1_EO)); % Create and empty array to put the filtered time series

        % Obtain frequency peaks
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

        %% HOPF MODEL
        
        % Load the SC matrix ans set the alphas based on the diagnostic group
        SC_selected = load([dirDTI '/' order{nreg} '.mat']);
        name = [selected_group '_' order{nreg} '_structure'];
        SC = SC_selected.structure;
        SC = SC/max(max(SC))*0.2; % Scale SC matrix to 0.2 

        %load alpha parameters and add the delta to simlate each condition
        load([dirA  order{nreg} '.mat']);
        a_vector = mean(a_vector);
        a_vector_perturb = a_vector + delta_a;

        % Parameters for the Hopf Model
        bp = 0;
        N = length(SC);

        disp('############processing#############')
        disp([order{nreg} selected_group drugSim(drug)])
        disp(datetime)

        SC = SC/max(max(SC))*0.2; % Scale SC matrix to 0.2 
        IntegStims = NaN(299,nSim,90); 
        TR = 1/fs; % Repetition time (1/fs)
        omega = repmat(2*pi*f_peaks,1,2);  % Changes the frequencies from Hz to radians and copies the array to a second dimension
        omega(:,1) = -omega(:,1); % Change the the sign of the first column to the opposite 
        dt=0.1.*(TR/2); % Compute the time steps from the repetition time to generate the simulated data
        sigma = 0.04; % Fixed value of sigma. In this paper it's 0.04
        dsig = sigma*sqrt(dt); % Scale sigma to the size of the time steps in the simulation?
        Nnew = size(SC,1); % Number of regions (ROIs) in the SC matrix. 90 in the AAL atlas
        isd = find(tril(ones(Nnew),-1)); % This tells us the indexes of the matrix that are in the lower diagonal
        StimStepTR=400; % Time to be simulated 
        StimStep=StimStepTR*TR;
        tmax=StimStepTR*nSim; 
        dt=0.1*TR/2;
        Tspan5000s=0:dt:((tmax-1)*TR);
        Tspan1Sim=0:dt:((StimStepTR-1)*TR);
        %allocte resting state netowrk alpha indexes
        grouping_prior= grouping.data_ICA_9;

        % No perturbation
        [IntegStim, integBase, meanIB, maxIB, minIB] = noperturbsim_vector_v2_non_2TR(G_perturb,SC,bp,StimStep,StimStepTR,Tspan1Sim,omega,TR,nSim,grouping_prior,a_vector_perturb);
        % Iterate over all ROIS (90)
        parfor pInd = 1:size(SC,1)

            % Synchronization perturbation
           [IntegStims(:,:,pInd), evIntegs(:,pInd), errIntegs(:,pInd), PILI(:,pInd), PILI_norm(:,pInd)]=perturbsynchsim_PILI_vectorized_non_2TR(G_perturb,SC,bp,StimStep,StimStepTR,Tspan1Sim,omega,pInd,TR,integBase,nSim,grouping_prior,a_vector_perturb);

        end

        %% Figures
        % Synch protocol
        options.x_axis = 0:1:StimStepTR-100-2; 
        options.handle     = figure(nFig);
        options.color_area = [128 193 219]./255;    
        options.color_line = [ 52 148 186]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'std';

        plot_areaerrorbar(squeeze(mean(IntegStims,2)).',options)

        hold on

        nFig = nFig + 1;

        % Draw line with basal state integration
        green = [0.4660 0.6740 0.1880];
        IBline = meanIB.*ones(1,size(IntegStim,1));
        plot(options.x_axis,IBline,'Linewidth',options.line_width,'color',green);
        xlim([0 300])
        title([order{nreg} ' LSD ' drugSim{drug} ' Synch'])
        ylabel('Integration')
        xlabel('Time (seconds)')
        set(gca,'fontsize', 14)

        mkdir(dirSaveFigures_synch)
        savefig([dirSaveFigures_synch order{nreg} ' LSD '  drugSim{drug} '.fig'])


        % Save the perturbation results
        integration_base = IntegStim;
        integration_base_mean = meanIB;

        integration_synch = squeeze(mean(IntegStims,2)).';

        mkdir(dirSave)
        
        save([dirSave order{nreg} ' ' selected_group ' ' drugSim{drug}], 'integration_base', 'integration_base_mean', 'integration_synch', 'PILI', 'PILI_norm')

        clearvars G optimalG IntegStim IntegStimn IntegStims W1_EO FC f_peaks_ W1_EO_corr evIntegs evIntegn errIntegs errIntegn

    end
end