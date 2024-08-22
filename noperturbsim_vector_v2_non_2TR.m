function [IntegStim, integbase, meanintegbase, maxintegbase, minintegbase]=noperturbsim_vector_v2(wG,C,bp,StimStep,StimStepTR,tsp,omega,TR, nsim,  grouping, a_vector)

dt=0.1*TR/2;
tmax=50; 
Tspan3000s=0:dt:tmax;
sig=0.04;
dsig = sqrt(dt)*sig; 
N=length(C);

%FILT


delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = .04;           % lowpass frequency of filter
fhi = .07;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter

for sim_n = 1:nsim
    
    xs=zeros(tmax,N);
    
    wC = wG*C;
    sumC = repmat(sum(wC,2),1,2);% for sum Cij*xj
    a_original = [grouping*a_vector.' grouping*a_vector.'];
    a = a_original;
    
    z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    
    % Get the time points where z should be sampled
    t_samp = tsp;
    sampling_vec = 1:length(t_samp);
    interval = sampling_vec(abs(mod(t_samp,TR))<0.01); % These are the points in t where the condition is met
    interval = interval(2)-interval(1); % Get the first two,samples where this happens, now we now the sampling interval
    
    % Warm up the Hopf model
    for t=0:dt:3000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    
    count = interval; % Sample at the first sample
    
    for t=tsp

        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        
        if count == interval % If it's time to sample
            nn = nn+1;
            xs(nn,:)=z(:,1)'; % The component that we want is the first one (real part of the complex number, which is the one that resembles fMRI dynamics)
            count = 0; % Reset the sampling counter
        end
        
        count = count + 1;
        
    end
    
    %%%%%%%%%%%%
    
    for seed=1:N
        ts=detrend(demean(xs(:,seed)'));
        tss = filtfilt(bfilt2,afilt2,ts);
        Xanalytic = hilbert(demean(tss));
        Phases(seed,:) = angle(Xanalytic);
    end
    
    T=1:1:size(Phases,2);
    integb=zeros(1,length(T));
    for t=T
        phasematrix=zeros(N,N);
        for i=1:N
            for j=1:N
                phasematrix(i,j)=exp(-3*adif(Phases(i,t),Phases(j,t)));
            end
        end
        cc=phasematrix;
        cc=cc-eye(N);
        pp=1;
        PR=0:0.01:0.99;
        cs=zeros(1,length(PR));
        for p=PR
            A=abs(cc)>p;
            [~, csize]=get_components(A);
            cs(pp)=max(csize);
            pp=pp+1;
        end
        integb(t)=sum(cs)*0.01/N;
    end
    
    % Assign values fort each simulation
    IntegStim_temp = integb;
    
    %%%%%%%%%%%%%%%%%%%%%
    
    initTime = 101; % Get the values after the first 100 seconds because that's when the perturbation is applied
    IntegStim(:,sim_n) = IntegStim_temp(initTime+1:end);
end

    integbase=mean(IntegStim,2);
    meanintegbase=mean(integbase);
    maxintegbase=max(integbase);
    minintegbase=min(integbase);

end

