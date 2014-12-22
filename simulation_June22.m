%clear all, example code
stream = RandStream('twister','Seed', 2014);
reset(stream);  % set up a seed for simulation
addpath(genpath(cat(2, pwd, '/PACErelease2.11')))

%% Generate simulation data

n = 100; % Number of curves
p = 40; % Numberof pooled grid points
t = (0 : (pi/2)/(p-1) : (pi/2)); %Pooled grid
sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter
nu = 2.5;
dense = 0.5;

%% generate stationary gaussian data

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

%% Non-stationary Gaussian data

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

%% common grid
nu = 2.5;  Burnin = 2000; M = 10000;  delta = 30; 

b = 1;  bs = 20; cgrid = 2; mat = 2;
%bhmcmc_fix_rho(Y_cgrid, T_cgrid, nu, delta, cgrid, b, bs,  Burnin, M, mat)
BHM_cgrid = bhmcmc_fix_rho(Y_cgrid, T_cgrid, nu, delta, 2, b, bs,  Burnin, M, mat); 
BHM_sparse = bhmcmc_fix_rho(Y, T, nu, delta, 0, 10, 500, Burnin, M, 2); 

BHMns_cgrid = bhmcmc_fix_rho(Yns_cgrid, Tns_cgrid, nu, delta, cgrid, b, 10,  Burnin, M, mat); 
BHMns_sparse = bhmcmc_fix_rho(Yns, Tns, nu, delta, 0, 10, 500, Burnin, M, 2); 

   
param_X = setOptions('regular', cgrid,'FVE_threshold',1,'corrPlot',0,'rho',-1); 
PCA_cgrid= FPCA(Y_cgrid, T_cgrid, param_X);   %perform PCA on x, common grid case
PCAns_cgrid= FPCA(Yns_cgrid, Tns_cgrid, param_X);   %perform PCA on x, common grid case

param_X_sparse = setOptions('regular', 0, 'FVE_threshold', 0.99, 'corrPlot', 0, 'rho', -1); 
PCA_sparse= FPCA(Y, T, param_X_sparse);   %perform PCA on x, sparse case
PCAns_sparse= FPCA(Yns, Tns, param_X_sparse);   %perform PCA on x, sparse case

%%

BHMns_sparse = bhmcmc_fix_rho(Yns, Tns, nu, delta, 0, 10, 500, Burnin, M, 2); 

save('simulationJune22.mat', 'GausData', 'GausData_NS', 'BHM_cgrid', 'PCA_cgrid', ...
            'BHMns_cgrid', 'PCAns_cgrid', 'BHM_sparse', 'BHMns_sparse', ...
            'PCA_sparse', 'PCAns_sparse')




