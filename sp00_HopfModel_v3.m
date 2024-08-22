function [FC_sim, FC_sim_tau, ku] = sp00_HopfModel_v2(wC, sumC, a, omega, dsig, Nnew, low_f, high_f, TTT, TR)

dt=0.1.*(TR/2);
xs=zeros(TTT,Nnew);
z = 0.1*ones(Nnew,2); % --> x = z(:,1), y = z(:,2)
nn=0;

%Tspan=0:dt:TTT;

% Get the time points where z should be sampled
t_samp = 0:dt:((TTT)*TR);
sampling_vec = 1:length(t_samp);
interval = sampling_vec(abs(mod(t_samp,TR))<0.01); % These are the points in t where the condition is met
interval = interval(2)-interval(1); % Get the first two,samples where this happens, now we now the sampling interval
%sampling_points = 1:interval:sampling_vec(end);

% Warm up the Hopf model
for t=0:dt:3000
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(Nnew,2);
end

count = interval; % Sample at the first sample

for t=0:dt:(((TTT)*TR))-1
    
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(Nnew,2);
    %z_complete(count,:) = z(:,1)';

    if count == interval % If it's time to sample
        nn = nn+1;
        xs(nn,:)=z(:,1)'; % The component that we want is the first one (real part of the complex number, which is the one that resembles fMRI dynamics, see line 5)
        count = 0; % Reset the sampling counter
    end
    
    count = count + 1;
    
end

% Filter simulated data
[ts, ku] = sp03_BandpassSimulatedTS(xs', low_f, high_f, TR);

FC_sim = corrcoef(ts');
t1 = ts(:,1:end-1);
t2 = ts(:,2:end);
FC_sim_tau = 1 - pdist2(t1,t2,'correlation');

% Compute FCD - very slow

for seed = 1:Nnew
    %     H(seed,:) = hilbert(ts(seed,:));
    H(seed,:) = hilbert(sp04_demean(ts(seed,:)));
end
ku = abs(sum(complex(cos(angle(H)),sin(angle(H))))/Nnew);
ku = ku(11:end-10);

end