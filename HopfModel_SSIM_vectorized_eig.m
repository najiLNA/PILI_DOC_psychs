function [euc_distance] = HopfModel_SSIM_vectorized_eig(wC, sumC, a_vector, omega, dsig, Nnew, low_f, high_f, TTT, TR, grouping, FEmp)

dt=0.1.*(TR/2);
xs=zeros(TTT,Nnew);
z = 0.1*ones(Nnew,2); % --> x = z(:,1), y = z(:,2)
nn=0;

Tspan=0:dt:TTT;

% Set the a paramaters according to the grouping into anatomical priors

a_1 = ((grouping)*a_vector.');

parfor population = 1:size(a_1,2)
    a = [a_1(:,population) a_1(:,population)];
    euc_distance(population) = pocket_hopf_eig(wC,z,sumC,a,dt,omega,Nnew,TR,TTT,dsig,nn,FEmp,low_f,high_f);
end

end