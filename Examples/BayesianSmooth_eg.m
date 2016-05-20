%% Set up seed for simulation 

stream = RandStream('twister','Seed', 2016);
reset(stream);  % set up a seed for simulation

%% add pathes of all source files of the required existing MATLAB libraries
    addpath(genpath(cat(2, pwd, '/bspline')))
    addpath(genpath(cat(2, pwd, '/fdaM'))) 
    addpath(genpath(cat(2, pwd, '/mcmcdiag'))) 
    addpath(genpath(cat(2, pwd, '/PACErelease2.11')))  

% load('./Examples/Data/Simu_Data.mat')
% load('./Examples/Data/Simu_Output.mat')

%% Set up parameters for simulation
n = 30; % Number of functional curves
p = 40; % Number of pooled grid points, or evaluated grid points
sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter in matern function
nu = 3.5; % order in matern function
pgrid = (0 : (pi/2)/(p-1) : (pi/2)); % Pooled grid
dense = 0.6; % proportion of observed data on the pooled grid
au = 0; bu = pi/2; % domain of t
m = 20; % Number of working grid points
stat = 1; cgrid = 1;

%% Stationary functional data on common grid
% generate simulated data
GausFD_cgrid = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, cgrid, stat);

% setup parameters for BFDA
% run with BHM
param = setOptions_bfda('smethod', 'bhm', 'cgrid', 1, 'mat', 1, ...
           'M', 10000, 'Burnin', 2000, 'w', 1, 'ws', 1);

% run with Bayesian Functional PCA
% param = setOptions_bfda('smethod', 'bfpca', 'M', 50, 'Burnin', 20);

% run with standard Bayesian Gaussion Process model
% param = setOptions_bfda('smethod', 'bgp', 'mat', 1, ...
%           'M', 50, 'Burnin', 20);

% run with Cubic Smoothing Splines
% param = setOptions_bfda('smethod', 'css', 'mat', 1, 'M', ...
%            50, 'Burnin', 20, 'pace', 0);

% call BFDA
[out_cgrid, param ] = ...
    BFDA(GausFD_cgrid.Xraw_cell, GausFD_cgrid.Tcell, param);

figure()
plot(out_cgrid.Z)
figure()
plot(out_cgrid.Sigma)

%% Stationary functional data on uncommon grid
GausFD_ucgrid = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, 0, stat);

param_uc = setOptions_bfda('smethod', 'bhm', 'cgrid', 0, 'mat', 1, 'M',...
    10000, 'Burnin', 2000, 'pace', 1, 'ws', 0.1);

[out_ucgrid, param_uc] = ...
    BFDA(GausFD_ucgrid.Xraw_cell, GausFD_ucgrid.Tcell, param_uc);

figure()
plot(out_ucgrid.Z)
figure()
plot(out_ucgrid.Sigma)

%% Non-stationary functional data on common grid
GausFD_cgrid_ns = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, cgrid, 0);

param_ns = setOptions_bfda('smethod', 'bhm', 'cgrid', 1, 'mat', 0, 'M',...
    10000, 'Burnin', 2000, 'pace', 1, 'ws', 0.01);

[out_cgrid_ns, param_ns] = ...
    BFDA(GausFD_cgrid_ns.Xraw_cell, GausFD_cgrid_ns.Tcell, param_ns);

figure()
plot(out_cgrid_ns.Z)
figure()
plot(out_cgrid_ns.Sigma)

%% Non-stationary functional data on uncommon grid
GausFD_ucgrid_ns = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, 0, 0);

param_uc_ns = setOptions_bfda('smethod', 'bhm', 'cgrid', 0, 'mat', 0, ...
    'M', 10000, 'Burnin', 2000, 'pace', 1, 'ws', 0.01);

[out_ucgrid_ns, param_uc_ns ] = ...
    BFDA(GausFD_ucgrid_ns.Xraw_cell, GausFD_ucgrid_ns.Tcell, param_uc_ns);

figure()
plot(out_ucgrid_ns.Z)
figure()
plot(out_ucgrid_ns.Sigma)

%% Stationary functional data on random grids
GausFD_rgrid = sim_gfd_rgrid(n, m, au, bu, sf, snr, nu, rho, stat);

param_rgrid = setOptions_bfda('smethod', 'babf', 'cgrid', 0, 'mat', 1, ...
    'M', 10000, 'Burnin', 2000, 'm', m, 'eval_grid', pgrid, 'ws', 1, ...
    'trange', [au, bu]);

% call BFDA
[out_rgrid, param_rgrid]= ...
    BFDA(GausFD_rgrid.Xraw_cell, GausFD_rgrid.Tcell, param_rgrid);

figure()
plot(out_rgrid.Z_cgrid)
figure()
plot(out_rgrid.Sigma_cgrid)

%% Nonstationary functional data on random grids
GausFD_rgrid_ns = sim_gfd_rgrid(n, m, au, bu, sf, snr, nu, rho, 0);

param_rgrid_ns = setOptions_bfda('smethod', 'babf', 'cgrid', 0, 'mat', ...
    0, 'M', 10000, 'Burnin', 2000, 'm', m, 'eval_grid', pgrid, 'ws', 0.05, ...
    'trange', [au, bu]);

% call BFDA
[out_rgrid_ns, param_rgrid_ns] = ...
    BFDA(GausFD_rgrid_ns.Xraw_cell, GausFD_rgrid_ns.Tcell, param_rgrid_ns);

figure()
plot(out_rgrid_ns.Z_cgrid)
figure()
plot(out_rgrid_ns.Sigma_cgrid)

%% Calculate rmse 
rmse(out_cgrid.Sigma_SE, GausFD_cgrid.Cov_true)
Xtrue_mat = reshape(cell2mat(GausFD_cgrid.Xtrue_cell), p, n);
rmse(out_cgrid.Z, Xtrue_mat)

% calculate the true non-stationary covariance matrix
Ctrue_ns = cov_ns(pgrid, sf, nu, rho);
rmse(out_cgrid_ns.Sigma_SE, Ctrue_ns)
    
%% Save simulated data sets
save('./Examples/Data/Simu_Data.mat', 'GausFD_cgrid', 'GausFD_ucgrid', ...
                           'GausFD_cgrid_ns', 'GausFD_ucgrid_ns', ...
                           'GausFD_rgrid', 'GausFD_rgrid_ns')

save('./Examples/Data/Simu_Output.mat', 'out_cgrid', 'out_ucgrid', ...
                           'out_cgrid_ns', 'out_ucgrid_ns', ...
                           'out_rgrid', 'out_rgrid_ns')
                       
%% Make plots by calling plot_script_smooth.m
close all;
run('./Examples/plot_script_smooth.m');

