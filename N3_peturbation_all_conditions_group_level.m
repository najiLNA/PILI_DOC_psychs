% This script uses the optimal G and a values to build models fit at the group level to
% defined strucutral and functional connectivities for each condition and then calculates PILI
% through perturbing each node of the AAL. 

clear
close all


fc1 = 0.04; % Lower bound for the band-pass filter
fc2 = 0.07; % Upper bound for the band-pass filter

nSim = 100; % Number of times to simulate the perturbation
IntegStims = NaN(299,nSim,90); 
nFig = 1;

dirG = [pwd '\Optimal G states\'];
dirA = [pwd '\Optimal Alphas states\'];

dirSave = [pwd '\perturbation group states'];
dirSaveFigures =[pwd '\perturbation group states'];

dirDTI = [pwd '\structural_data\'];
dirDTI_MCS = [dirDTI '\MCS_structure_mat_files\'];
dirDTI_UWS = [dirDTI '\UWS_structure_mat_files\'];
dirDTI_CNT = [dirDTI '\avg\'];

recordings_CNT=dir([dirDTI_CNT '*.mat']);
recordings_MCS=dir([dirDTI_MCS '*.mat']);
recordings_UWS=dir([dirDTI_UWS '*.mat']);

drugSim = {'LSD','LSD_PCB', 'PSILO_maas', 'PSILO_maas_PCB', 'PSILO_imp', 'PSILO_imp_PCB', 'Propofol_W', 'Propofol_S1', 'Propofol_S2' , 'MCS', 'UWS', 'CNT', 'Dex_W', 'Dex_S1', 'Dex_S2'};

grouping = load([ pwd '\data_ICA_9.mat']);
grouping.labels = {'ATTENTION', 'AUD', 'DMN', 'FP', 'HPC', 'PRECUNEUS', 'SM', 'THAL','VIS'};


for drug=1:length(drugSim)
    drugSim{drug}
    % Loads the parameters for each condition 
    switch drugSim{drug}
            case 'LSD'
                dataDir_key = 'Data_LSD\LSD\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\LSD\LSD_Rest1_denoise2.mat']);
                a_vector = load([dirA '\LSD\LSD_Rest1_denoise2.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
                
            case 'LSD_PCB'
                dataDir_key = 'Data_LSD\PCB\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\LSD\PCB_Rest1_denoise2.mat']);
                a_vector = load([dirA '\LSD\PCB_Rest1_denoise2.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;

            case 'PSILO_maas'
                dataDir_key = 'Data_PSILO\maas\PSILO\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\psilo\maas\PSILO\data_PSILO_maas.mat']);
                a_vector = load([dirA '\PSILO_maas\data_PSILO_maas.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/1.3;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
            case 'PSILO_maas_PCB'
                dataDir_key = 'Data_PSILO\maas\PCB\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\psilo\maas\data_PCB_maas.mat']);
                a_vector = load([dirA '\PSILO_maas\data_PSILO_maas.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/1.3;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
            case 'PSILO_imp'
                dataDir_key = 'Data_PSILO\imperial\PSILO\';
                G = load([dirG '\psilo\imperial\PSILO_2nd_half_ica.mat']);
                a_vector = load([dirA '\PSILO_imperial\PSILO_2nd_half_ica.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/3;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
            case 'PSILO_imp_PCB'
                dataDir_key = 'Data_PSILO\imperial\PCB\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\psilo\imperial\PCB_2nd_half_ica.mat']);
                a_vector = load([dirA '\PSILO_imperial\PCB_2nd_half_ica.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/3;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
             case 'Propofol_W1' 
                dataDir_key = 'Data_propofol\W1\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\propofol\W1_.mat']);
                a_vector = load([dirA '\propofol\W1_.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2.4;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
                
             case 'Propofol_S1' 
                dataDir_key = 'Data_propofol\S1\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\propofol\S1_.mat']);
                a_vector = load([dirA '\propofol\S1_.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2.4;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
                
             case 'Propofol_S2' 
                dataDir_key = 'Data_propofol\S2\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir '*.mat']);
                G = load([dirG '\propofol\S2_.mat']);
                a_vector = load([dirA '\propofol\S2_.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2.4;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
             case 'MCS' 
                dataDir_key = 'Data_DoC\';
                dataDir = [pwd '\BOLD_data\' dataDir_key 'MCS'];
                recordings=dir([dataDir 'DoC_MCS.mat']);
                G = load([ dirG '\Optimal G DoC groups\DoC_MCS.mat']);
                a_vector = load([dirA '\DoC\DoC\MCS.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2;
                load([dirDTI_MCS  'MCS_avg_structure.mat']);
                SC=MCS_avg_structure;
              case 'UWS' 
                dataDir_key = 'Data_DoC\';
                dataDir = [pwd '\BOLD_data\' dataDir_key '\UWS'];
                recordings=dir([dataDir 'DoC_UWS.mat']);
                G = load([ dirG '\Optimal G DoC groups\DoC_UWS.mat']);
                a_vector = load([dirA '\DoC\DoC\UWS.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2;
                load([dirDTI_UWS  'UWS_avg_structure.mat']);
                SC=UWS_avg_structure;
              case 'CNT' 
                dataDir_key = 'Data_DoC\';
                dataDir = [pwd '\BOLD_data\' dataDir_key];
                recordings=dir([dataDir 'DoC_CNT.mat']);
                G = load([ dirG '\Optimal G DoC groups\DoC_CNT.mat']);
                a_vector = load([dirA '\DoC\DoC\CNT.mat']);
                a_vector = mean(a_vector.a_vector);
                fs=1/2;
                load([dirDTI_CNT  'CNT_avg_structure.mat']);
                SC=CNT_avg_structure;
                
            
    end
    


    
        
    
        G=G.optimalG;
        genericMatrix = struct2cell(load([dataDir '\' recordings(1).name]));
        W1.time_series = genericMatrix{1};
       
        W1.name = drugSim{drug};

      
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
            % Obtain peak frequencies from the time series 
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
     

        % Parameters for the Hopf Model
        bp = 0;
        N = length(SC);

        disp('############processing#############')
        disp([drugSim(drug)])
        disp(datetime)


        SC = SC/max(max(SC))*0.2; % Scale SC matrix to 0.2 anyway

        TR = 1/fs; % Repetition time (1/fs)
        %Tmax = 300; % Time to be simulated in volumes. 300 for now, might change it to less time? Now it's 217 like the actual data
        omega = repmat(2*pi*f_peaks,1,2);  % Changes the frequencies from Hz to radians and copies the array to a second dimension
        omega(:,1) = -omega(:,1); % Change the the sign of the first column to the opposite (I guess this is needed for the vectorization of the model they have going on)
        dt=0.1.*(TR/2); % Compute the time steps from the repetition time to generate the simulated data

        sigma = 0.04; % Fixed value of sigma. In this paper it's 0.04
        dsig = sigma*sqrt(dt); % Scale sigma to the size of the time steps in the simulation?
        Nnew = size(SC,1); % Number of regions (ROIs) in the SC matrix. 90 in the AAL atlas
        isd = find(tril(ones(Nnew),-1)); % This tells us the indexes of the matrix that are in the lower diagonal (which has all the info since the SC matrix is symmetrical
        StimStepTR=400; % Time to be simulated in seconds. 
        StimStep=StimStepTR*TR;
        tmax=StimStepTR*nSim; 
        dt=0.1*TR/2;
        Tspan5000s=0:dt:((tmax-1)*TR);
        Tspan1Sim=0:dt:((StimStepTR-1)*TR);
        grouping_prior= grouping.data_ICA_9;

        % No perturbation
        [IntegStim, integBase, meanIB, maxIB, minIB] = noperturbsim_vector_v2_non_2TR(G,SC,bp,StimStep,StimStepTR,Tspan1Sim,omega,TR,nSim,grouping_prior,a_vector);
        % Iterate over all ROIS (90)

        parfor pInd = 1:size(SC,1)

            % perturbation
           [IntegStims(:,:,pInd), evIntegs(:,pInd), errIntegs(:,pInd), PILI(:,pInd), PILI_norm(:,pInd)]=perturbsynchsim_PILI_vectorized_non_2TR(G,SC,bp,StimStep,StimStepTR,Tspan1Sim,omega,pInd,TR,integBase,nSim,grouping_prior,a_vector);

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
        title([drugSim{drug} ' Synch'])
        ylabel('Integration')
        xlabel('Time (seconds)')
        set(gca,'fontsize', 14)

        mkdir(dirSaveFigures)
        savefig([dirSaveFigures drugSim{drug} '_Synch.fig'])


        % Save the perturbation results
        integration_base = IntegStim;
        integration_base_mean = meanIB;

        integration_synch = squeeze(mean(IntegStims,2)).';

        mkdir(dirSave)
        
        save([dirSave drugSim{drug}], 'integration_base', 'integration_base_mean', 'integration_synch', 'PILI', 'PILI_norm')

        clearvars G optimalG IntegStim IntegStimn IntegStims W1_EO FC f_peaks_ W1_EO_corr evIntegs evIntegn errIntegs errIntegn

    
end