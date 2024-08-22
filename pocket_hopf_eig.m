function euc_distance = pocket_hopf_eig(wC,z,sumC,a,dt,omega,Nnew,TR,TTT,dsig,nn,FEmp,low_f,high_f)

for t=0:dt:(((TTT)*TR)+3000.1)
    suma = wC*z - sumC.*z;
    zz = z(:,end:-1:1);
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(Nnew,2);
    
    if (abs(mod(t,TR))<0.01) && (t>3000) % If the temporal sample t corresponds to a TR, update the generated time series
        nn=nn+1;
        xs(nn,:)=z(:,1)'; % The component that we want is the first one (real part of the complex number, which is the one that resembles fMRI dynamics, see line 5)
    end
end

% Filter simulated data
[ts, ku] = sp03_BandpassSimulatedTS(xs', low_f, high_f, TR);

FSim = corrcoef(ts');
t1 = ts(:,1:end-1);
t2 = ts(:,2:end);
% Get lower diagonal
m  = (1:size(FSim,1)).' >= (1:size(FSim,2));
v_fsim  = FSim(m);
v_femp  = FEmp(m);

euc_distance = pdist2(v_fsim.',v_femp.','euclidean');

