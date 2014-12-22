%% Set up seed for simulation 
%clear all, example code
stream = RandStream('twister','Seed', 2014);
reset(stream);  % set up a seed for simulation

%% Set up parameters for simulation

n = 50; % Number of curves
p = 80; % Numberof pooled grid points
t = (0 : (pi/2)/(p-1) : (pi/2)); %Pooled grid
sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter
nu = 3.5;
dense = 0.6;

%% Generate stationary gaussian data

load('SimulatedFunctionalData.mat');
%GausData = gaussian_sim(t, n, sf, snr, nu, rho, dense);

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

%% Generate Non-stationary Gaussian data

%GausData_NS = gaussian_simNS(t, n, sf, snr, nu, rho, dense);

tns = GausData_NS.t;
Cns = GausData_NS.C;
Dns = GausData_NS.D;
X0ns = GausData_NS.X0;
Y0ns = GausData_NS.Y0;
mutns = GausData_NS.mut;


Xns = GausData_NS.X;
Yns = GausData_NS.Y;
Tns = GausData_NS.T;

Xns_cgrid = GausData_NS.X_cgrid;
Yns_cgrid = GausData_NS.Y_cgrid;
Tns_cgrid = GausData_NS.T_cgrid;

%% Run PACE with simulated data

addpath(genpath(cat(2, pwd, '/PACErelease2.11')))

% param_X = setOptions('regular', 1,'FVE_threshold',1,'corrPlot',0,'rho',-1); 
% PCA_cgrid= FPCA(Y_cgrid, T_cgrid, param_X);   %perform PCA on x, common grid case
    Phihat=getVal(PCA_cgrid,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_cgrid,'lambda'); % The estimated eigenvalues.
    Sigma_est = Phihat*diag(lamhat)*Phihat';
    mu_cgrid = (getVal(PCA_cgrid, 'mu'));
    
% PCAns_cgrid= FPCA(Yns_cgrid, Tns_cgrid, param_X);   %perform PCA on x, common grid case
    Phihat_ns=getVal(PCAns_cgrid,'phi');  % The estimated Eigenfunctions.
    lamhat_ns=getVal(PCAns_cgrid,'lambda'); % The estimated eigenvalues.
    Sigma_est_ns = Phihat_ns * diag(lamhat_ns) * Phihat_ns';
    mu_ns = (getVal(PCAns_cgrid, 'mu'));
    
    
%param_X_sparse = setOptions('regular', 0, 'FVE_threshold', 0.99, 'corrPlot', 0, 'rho', -1); 
% PCA_sparse= FPCA(Y, T, param_X_sparse);   %perform PCA on x, sparse case
    Phihat_sparse=getVal(PCA_sparse,'phi');  % The estimated Eigenfunctions.
    lamhat_sparse=getVal(PCA_sparse,'lambda'); % The estimated eigenvalues.
    Sigma_est_sparse = Phihat_sparse * diag(lamhat_sparse) * Phihat_sparse';
    mu_sparse = (getVal(PCA_sparse, 'mu'));
    
    
% PCAns_sparse= FPCA(Yns, Tns, param_X_sparse);   %perform PCA on x, sparse case
    Phihat_sparse_ns=getVal(PCAns_sparse,'phi');  % The estimated Eigenfunctions.
    lamhat_sparse_ns=getVal(PCAns_sparse,'lambda'); % The estimated eigenvalues.
    Sigma_est_sparse_ns = Phihat_sparse_ns * diag(lamhat_sparse_ns) * Phihat_sparse_ns';
    mu_sparse_ns = (getVal(PCAns_sparse, 'mu'));

%% common grid
Burnin = 2000; M = 10000;  delta = 5; 
b = 20; mat = true;
 
%%
    BHM_cgrid = bhmcmc(Y_cgrid, T_cgrid, delta, 2, ...
        Burnin, M, mat, Sigma_est, mu_cgrid, [], 7.2); 
    rmse(BHM_cgrid.iK, C) 
    rmse(Sigma_est, C)
    rmse(COR(BHM_cgrid.iK), COR(C))
    rmse(COR(Sigma_est), COR(C))
    rmse(BHM_cgrid.Z, X0)
    
    BHMns_cgrid = bhmcmc(Yns_cgrid, Tns_cgrid, delta, 2, ...
        Burnin, M, 0, Sigma_est_ns, mu_ns, [], 5); 
    rmse(BHMns_cgrid.iK, Cns)
    rmse(Sigma_est_ns, Cns)
    rmse(COR(BHMns_cgrid.iK), COR(Cns))
    rmse(COR(Sigma_est_ns), COR(Cns))
    rmse(BHMns_cgrid.Z, X0ns)

    
        BHM_sparse = bhmcmc(Y, T, delta, 0, Burnin, M, mat, ...
        Sigma_est_sparse, mu_sparse, [], 1); 
    rmse(BHM_sparse.iK, C)
    rmse(Sigma_est_sparse, C)
    rmse(BHM_sparse.Z, X0)
    
    BHMns_sparse = bhmcmc(Yns, Tns, delta, 0, Burnin, M, 0,...
        Sigma_est_sparse_ns, mu_sparse_ns, [], 1); 
    rmse(BHMns_sparse.iK, Cns)
    rmse(Sigma_est_sparse_ns, Cns)
        rmse(BHMns_sparse.Z, X0ns)

    

    
% rmse(BHMns_cgrid.iK, Cns)
   
% bhmcmc(Y, T, delta, cgrid, b, bs, Burnin, M, mat, Sigma_est, mu_p, nu)

%%
%July5

save('simulationAug8.mat', 'GausData', 'GausData_NS', 'BHM_cgrid', 'PCA_cgrid', ...
            'BHMns_cgrid', 'PCAns_cgrid', 'BHM_sparse', 'BHMns_sparse', ...
            'PCA_sparse', 'PCAns_sparse')
        
 save('temp_cgrid.mat', 'BHM_cgrid', 'BHMns_cgrid')
 save('temp_sparse.mat', 'BHM_sparse', 'BHMns_sparse')
 



