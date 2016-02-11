%% Set up seed for simulation 

stream = RandStream('twister','Seed', 2016);
reset(stream);  % set up a seed for simulation

%% Set up parameters for simulation

n = 10; % Number of functional curves
p = 30; % Number of pooled grid points
pgrid = (0 : (pi/2)/(p-1) : (pi/2)); % Pooled grid
sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter in matern function
nu = 3.5; % order in matern function
dense = 0.6; % proportion of observed data on the pooled grid

load('./Examples/Data/SimulationData.mat')

%% Stationary functional data on common grid
% generate simulated data
GausFD_cgrid = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, cgrid, stat);

% setup parameters for BFDA
param = setOptions_bfda('cgrid', 1, 'mat', 1, 'babf', 0, ...
    'M', 50, 'Burnin', 20);

% call BFDA
out_cgrid =  BFDA(GausFD_cgrid.Y, GausFD_cgrid.T, param);


%% Stationary functional data on uncommon grid
GausFD_ucgrid = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, 0, stat);

param_uc = setOptions_bfda('cgrid', 0, 'mat', 1, 'babf', 0, ...
    'M', 50, 'Burnin', 20, 'pace', 0);

out_ucgrid =  BFDA(GausFD_ucgrid.Y, GausFD_ucgrid.T, param_uc);


%% Non-stationary functional data on common grid
GausFD_cgrid_ns = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, cgrid, 0);

param_ns = setOptions_bfda('cgrid', 1, 'mat', 0, 'babf', 0, ...
    'M', 50, 'Burnin', 20, 'pace', 1);

out_cgrid_ns =  BFDA(GausFD_cgrid_ns.Y, GausFD_cgrid_ns.T, param_ns);


%% Non-stationary functional data on uncommon grid
GausFD_ucgrid_ns = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, 0, 0);

param_uc_ns = setOptions_bfda('cgrid', 0, 'mat', 0, 'babf', 0, ...
    'M', 50, 'Burnin', 20, 'pace', 1);

out_ucgrid_ns =  BFDA(GausFD_ucgrid_ns.Y, GausFD_ucgrid_ns.T, param_uc_ns);


%% Calculate rmse 

    rmse(out_cgrid.iKSE, GausFD_cgrid.C)
    rmse(out_cgrid.Z, GausFD_cgrid.X0)
    

%% Save simulated data sets

save('./Examples/Data/SimulationData.mat', 'GausFD_cgrid', 'GausFD_ucgrid', ...
                           'GausFD_cgrid_ns', 'GausFD_ucgrid_ns')


