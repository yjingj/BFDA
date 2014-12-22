% cd /Users/jingjingyang/Dropbox/JingjingStuff/BayesSmoothing/BS_Matlab_Code
clear all, 
stream = RandStream('twister','Seed', 2014);
reset(stream);  % set up a seed for simulation

%% Generate simulation data

n = 100; % Number of curves
p = 40; % Numberof pooled grid points
t = (0 : (pi/2)/(p-1) : (pi/2)); %Pooled grid
sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter
nu = 2.5;
dense = 0.4;

%%
% GausData = gaussian_sim(t, n, sf, snr, nu, rho, dense);

t = GausData.t;
C = GausData.C;
D = GausData.D;
X0 = GausData.X0;
Y0 = GausData.Y0;
mut = GausData.mut;


X = GausData.X;
Y = GausData.Y;
T = GausData.T;

X_cgrid = GausData.X_cgrid;
Y_cgrid = GausData.Y_cgrid;
T_cgrid = GausData.T_cgrid;

% Y = Y_cgrid; T = T_cgrid; 


%% common grid
nu = 2.5;  Burnin = 2000; M = 10000; 
b = 1; cgrid = 2; delta = 30; mat = 2;

bs = 5;
BHM_cgrid = bhmcmc_fix_rho(Y_cgrid, T_cgrid, nu, delta, cgrid, b, bs,  Burnin, M, mat); 

% wrho = 0.03;
% Matern_cgrid = maternmod_mcmc(Y_cgrid, T_cgrid, nu, cgrid, b, bs, wrho, Burnin, M);

% Tprior_cgrid = tprior_mcmc(Y_cgrid, T_cgrid, delta, cgrid, b, Burnin, M);

%wrho_prop = 0.08;
% BHM_cgrid = bhmcmc(Y_cgrid, T_cgrid, nu, delta, cgrid, b, bs, wrho_prop, Burnin, M); 

addpath(genpath(cat(2, pwd, '/mcmcdiag')))


%%
cgrid = 2;
addpath(genpath(cat(2, pwd, '/PACErelease2.11')))
param_X = setOptions('regular', cgrid,'FVE_threshold',1,'corrPlot',0,'rho',-1); 
PCA_cgrid= FPCA(Y_cgrid, T_cgrid, param_X);   %perform PCA on x, common grid case

Phihat=getVal(PCA_cgrid,'phi');  % The estimated Eigenfunctions.
snhat2=getVal(PCA_cgrid,'sigma'); % The estimated noise for X(t)
lamhat=getVal(PCA_cgrid,'lambda'); % The estimated eigenvalues.
pcx = getVal(PCA_cgrid,'xi_est');  % The estimated PC scores. 

PCA_mu = getVal(PCA_cgrid, 'mu');
PCA_Z = pcx * Phihat';
Sigma_pace_cgrid=Phihat*diag(lamhat)*Phihat'; % The estimated covariance


%%  Comparing our Bayesian method

rimse(t, X0, BHM_cgrid.Z, 0)
% rimse(t, X0, Tprior_cgrid.Z, 0)
% rimse(t, X0, Matern_cgrid.Z, 0)
rimse(t, X0, PCA_Z' + repmat(PCA_mu', 1, 100), 0)

rimse(t, mut, BHM_cgrid.mu', 0)
% rimse(t, mut, Tprior_cgrid.mu', 0)
% rimse(t, mut, Matern_cgrid.mu', 0)
rimse(t, mut, PCA_mu, 0)

[xgrid,ygrid]=meshgrid(GausData.t); % 3D plots of the covariance matrices
figure()
mesh(xgrid,ygrid, Sigma_pace_cgrid)
figure()
mesh(xgrid, ygrid, BHM_cgrid.iK)
figure()
mesh(xgrid, ygrid,  C)

rimse(t, C, BHM_cgrid.iK, 2)
rimse(t, C, Tprior_cgrid.iK, 2)
rimse(t, C, Matern_cgrid.iK, 2)
rimse(t, C, Sigma_pace_cgrid, 2)

R = C ./ (sqrt(diag(C)) * sqrt(diag(C))' );
BHM_R = BHM_cgrid.iK ./ (sqrt(diag(BHM_cgrid.iK)) * sqrt(diag(BHM_cgrid.iK))' );
%Tprior_R = Tprior_cgrid.iK ./ (sqrt(diag(Tprior_cgrid.iK)) * sqrt(diag(Tprior_cgrid.iK))' );
%Matern_R = Matern_cgrid.iK ./ (sqrt(diag(Matern_cgrid.iK)) * sqrt(diag(Matern_cgrid.iK))' );
pace_R = Sigma_pace_cgrid ./ (sqrt(diag(Sigma_pace_cgrid)) * sqrt(diag(Sigma_pace_cgrid))' );

rimse(t, R, BHM_R, 2)
%rimse(t, R, Tprior_R, 2)
%rimse(t, R, Matern_R, 2)
rimse(t, R, pace_R, 2)


%%  Sparse Case
cgrid = 0;   
b = 1/3; bs = 50;
BHM_sparse = bhmcmc_fix_rho(Y, T, nu, delta, cgrid, b, bs, Burnin, M, mat); 
rimse(t, C, BHM_sparse.iK, 2)

param_X_sparse = setOptions('regular', cgrid, 'FVE_threshold', 0.99, 'corrPlot', 0, 'rho', -1); 
PCA_sparse= FPCA(Y, T, param_X_sparse);   %perform PCA on x, sparse case

pcx_sparse = getVal(PCA_sparse, 'xi_est');  % The estimated PC scores. 
Phihat_sparse=getVal(PCA_sparse,'phi');  % The estimated Eigenfunctions.
snhat2_sparse=getVal(PCA_sparse,'sigma'); % The estimated noise for X(t)
lamhat_sparse=getVal(PCA_sparse,'lambda'); % The estimated eigenvalues.

PCA_mu_sparse = getVal(PCA_sparse, 'mu');
PCA_Z_sparse = pcx_sparse * Phihat_sparse';
Sigma_pace_sparse=Phihat_sparse*diag(lamhat_sparse)*Phihat_sparse'; % The estimated covariance

%%  Comparing our Bayesian method

SCov = cov(Y0');
[xgrid,ygrid]=meshgrid(GausData.t);

figure()
mesh(xgrid, ygrid, SCov)
figure()
mesh(xgrid,ygrid, Sigma_pace_sparse)
figure()
mesh(xgrid, ygrid, BHM_sparse.iK)
figure()
mesh(xgrid, ygrid, C)

rimse(t, X0, BHM_sparse.Z, 0)
rimse(t, X0, PCA_Z_sparse', 0)

rimse(t, mut, BHM_sparse.mu', 0)
rimse(t, mut, PCA_mu_sparse, 0)

rimse(t, C, BHM_sparse.iK, 2)
rimse(t, C, Sigma_pace_sparse, 2)
rimse(t, C, SCov, 2)


%% Save outputs
save('StationaryGauss150p.mat', 'GausData', 'BHM_cgrid', 'PCA_cgrid', ...
           'BHM_sparse', 'PCA_sparse')
 
       
%save('StationaryGauss.mat', 'GausData', 'BHM_cgrid', 'PCA_cgrid', 'Matern_cgrid',...
  %         'BHM_sparse', 'PCA_sparse', 'Tprior_cgrid')
 






