%clear all;
stream = RandStream('twister','Seed', 2014);
reset(stream);  % set up a seed for simulation

%% Parameter set up
addpath(genpath(cat(2, pwd, '/PACErelease2.11')))
n = 50; % Number of curves
p = 80; % Numberof pooled grid points
t = (0 : (pi/2)/(p-1) : (pi/2)); %Pooled grid

sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter
nu = 3.5;
dense = 0.6; % percentage of data density on pooled grid
% Repeat the simulation for 100 times

Burnin = 2000; M = 10000;
b = 1; delta = 5;  cgrid = 2;

param_X = setOptions('regular', 2,'FVE_threshold',0.99,'corrPlot',0,'rho',-1); 
param_X_sparse = setOptions('regular', 0, 'FVE_threshold', 0.99, 'corrPlot', 0, 'rho', -1); 


%%
Nsim = 100;
YSim = NaN(p, n, Nsim);
XSim = NaN(p, n, Nsim);

YSim_NS = NaN(p, n, Nsim);
XSim_NS = NaN(p, n, Nsim);

% save results form BHM
ZSim = NaN(p, n, Nsim); % Smoothed curve
iKSim = NaN(p, p, Nsim); % Covariance matrix
muSim = NaN(p, Nsim); % Mean vector
rnSim = NaN(1, Nsim); % Precision of error, gamma_n
rsSim = NaN(1, Nsim); % signal variance, sf
rhoSim = NaN(1, Nsim); % gamma_theta1
nuSim = NaN(1, Nsim);

% save results from BHM_sparse
ZSim_NS = NaN(p, n, Nsim); % Smoothed curve
iKSim_NS = NaN(p, p, Nsim); % Covariance matrix
muSim_NS = NaN(p, Nsim); % Mean vector
rnSim_NS = NaN(1, Nsim); % Precision of error, gamma_n
rsSim_NS = NaN(1, Nsim); % signal variance, sf
rhoSim_NS = NaN(1, Nsim); % gamma_theta1
nuSim_NS = NaN(1, Nsim);

% save results from PCA
ZpSim = NaN(p, n, Nsim); % Smoothed curve
iKpSim = NaN(p, p, Nsim); % Covariance matrix
mupSim = NaN(p, Nsim); % Mean vector
snpSim = NaN(1, Nsim); % noise variance

ZpSim_NS = NaN(p, n, Nsim); % Smoothed curve
iKpSim_NS = NaN(p, p, Nsim); % Covariance matrix
mupSim_NS = NaN(p, Nsim); % Mean vector
snpSim_NS = NaN(1, Nsim); % noise variance


%%  

for (iSim = 1 : Nsim)

    %************************** Generate stationary Gaussian data
    GausData = gaussian_sim(t, n, sf, snr, nu, rho, dense);
    Y_cgrid = GausData.Y_cgrid;
    T_cgrid = GausData.T_cgrid;
    
    % simulated true and noise functional signals for non-sparse case
    XSim(:, :, iSim) = GausData.X0;
    YSim(:, :, iSim) = GausData.Y0;
    
    % PACE for non-stationary case
    PCA_cgrid= FPCA(Y_cgrid, T_cgrid, param_X);   %perform PCA on x, common grid case
    Phihat=getVal(PCA_cgrid,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_cgrid,'lambda'); % The estimated eigenvalues.
    pcx = getVal(PCA_cgrid,'xi_est');  % The estimated PC scores. 

    mupSim(:, iSim) = (getVal(PCA_cgrid, 'mu'))';
    snpSim(iSim) = getVal(PCA_cgrid,'sigma'); 
    ZpSim(:, :, iSim) = (pcx * Phihat')' + repmat(mupSim(:, iSim), 1, n);
    iKpSim(:, :, iSim) = Phihat*diag(lamhat)*Phihat'; % The estimated covariance

    BHM_cgrid = bhmcmc(Y_cgrid, T_cgrid, delta, 2, ...
        Burnin, M, 2, iKpSim(:, :, iSim), getVal(PCA_cgrid, 'mu'), []);  % BHM model  
    ZSim(:, :, iSim) = BHM_cgrid.Z; % Smoothed curve
    iKSim(:, :, iSim) = BHM_cgrid.iK; % Covariance matrix
    muSim(:, iSim) = BHM_cgrid.mu; % Mean vector 
    rnSim(:, iSim) = BHM_cgrid.rn; % Precision of error, gamma_n
    rsSim(:, iSim) = BHM_cgrid.rs; % signal variance, sf
    rhoSim(:, iSim) = BHM_cgrid.rho; % gamma_theta1
    nuSim(:, iSim) = BHM_cgrid.nu;

     %************************ Genrate Non-Stationary Gaussian data
    GausData_NS = gaussian_simNS(t, n, sf, snr, nu, rho, dense);
    Y_NS = GausData_NS.Y_cgrid;
    T_NS = GausData_NS.T_cgrid;
    
    XSim_NS(:, :, iSim) = GausData_NS.X0;
    YSim_NS(:, :, iSim) = GausData_NS.Y0;
    
    % PACE 
    PCA_NS= FPCA(Y_NS, T_NS, param_X);   %perform PCA on x, common grid case
    Phihat=getVal(PCA_NS,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_NS,'lambda'); % The estimated eigenvalues.
    pcx = getVal(PCA_NS,'xi_est');  % The estimated PC scores. 

    mupSim_NS(:, iSim) = (getVal(PCA_NS, 'mu'))';
    snpSim_NS(iSim) = getVal(PCA_NS,'sigma'); 
    ZpSim_NS(:, :, iSim) = (pcx * Phihat')' + repmat(mupSim_NS(:, iSim), 1, n);
    iKpSim_NS(:, :, iSim) = Phihat*diag(lamhat)*Phihat'; % The estimated covariance

    %BHM_cgrid model
    BHM_NS = bhmcmc(Y_NS, T_NS, delta, 2, ...
        Burnin, M, 0, iKpSim_NS(:, :, iSim), getVal(PCA_NS, 'mu'), []); 
    
    ZSim_NS(:, :, iSim) = BHM_NS.Z; % Smoothed curve
    iKSim_NS(:, :, iSim) = BHM_NS.iK; % Covariance matrix
    muSim_NS(:, iSim) = BHM_NS.mu; % Mean vector 
    rnSim_NS(:, iSim) = BHM_NS.rn; % Precision of error, gamma_n
    rsSim_NS(:, iSim) = BHM_NS.rs; % signal variance, sf


end



C = GausData.C;
mut = GausData.mut;
D = GausData.D;

C_NS = GausData_NS.C;
mut_NS = GausData_NS.mut;

save('StatSimAug10.mat', 'XSim', 'YSim', 'XSim_NS', 'YSim_NS', ...
     'ZSim', 'iKSim', 'muSim', 'rnSim', 'rsSim', 'rhoSim', ...
    'ZSim_NS', 'iKSim_NS', 'muSim_NS', 'rnSim_NS', 'rsSim_NS', 'rhoSim_NS', ...
    'ZpSim', 'snpSim', 'iKpSim', 'mupSim', 'ZpSim_NS', 'snpSim_NS', 'iKpSim_NS', 'mupSim_NS',...
    'C', 'mut', 'D', 'C_NS', 'mut_NS');


%% load('Device1_WRL.mat') real spectroscopy data from el paso 2013
%wcal 291 by 462

load('Device1_WRL.mat')
emWL = 410:700;

x = 410:3:700;
Yspecmat = wcal(1:3:291 , :)';

n = 462; p = 97;
Yspecleft1 = wcal( (1:96) .* 3 +1, :)';
xleft1 = emWL( (1:96) .* 3 +1 );

Yspecleft2 = wcal( (1:96) .* 3 +2, :)';
xleft2 = emWL((1:96) .* 3 +2);

Y = mat2cell(Yspecmat, ones(462, 1))';
T = mat2cell(repmat(x, 462, 1), ones(462, 1))';

% PACE 
addpath(genpath(cat(2, pwd, '/PACErelease2.11'))) 
param_X = setOptions('regular', cgrid,'FVE_threshold',1,'corrPlot',0,'rho',-1); 

   % PCA_spec= FPCA(Y, T, param_X);   %perform PCA on x, common grid case
    Phihat=getVal(PCA_spec,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_spec,'lambda'); % The estimated eigenvalues.
    pcx = getVal(PCA_spec,'xi_est');  % The estimated PC scores. 

    mu_spec_pace = (getVal(PCA_spec, 'mu'))';
    sn_spec_pace = getVal(PCA_spec,'sigma'); 
    Z_spec_pace = (pcx * Phihat')' + repmat(mu_spec_pace, 1, 462);
    iK_spec_pace = Phihat*diag(lamhat)*Phihat'; % The estimated covariance

    
delta = 5; cgrid = 2; b = 1; bs = 1; Burnin = 2000; M = 10000;

spec_BHM = bhmcmc(Y, T, delta, cgrid, b, bs, Burnin, M, 0, iK_spec_pace, mu_spec_pace');


 save('specAug10', 'Y', 'T', 'x', 'Yspecmat', 'spec_BHM', 'Z_spec_pace', 'mu_spec_pace', ...
           'iK_spec_pace', 'sn_spec_pace',  'PCA_spec',...
           'emWL', 'Yspecleft1', 'xleft1', 'Yspecleft2', 'xleft2')
   
       
   %% Obese study
load('sleep.mat')
load('ObeseJuly5.mat')

Sigma_est_ob = BHM_obese.Sigma_est;
mu_p_ob = BHM_obese.mu_p;

%addpath(genpath(cat(2, pwd, '/mcmcdiag')))
t = 100:200;

obese_data = sleep_ee(sleep_obese == 1, t);
non_obese_data = sleep_ee(sleep_obese == 0, t);

T = mat2cell(repmat(t, 44, 1), ones(1,44), length(t))';
Y = mat2cell(obese_data, ones(1, 44), length(t))';
BHM_obese = bhmcmc(Y, T, 5, 2, 1/1000, 1, 2000, 10000, 0, Sigma_est_ob, mu_p_ob, []); %


%Tnon = mat2cell(repmat(t, 62, 1), ones(1,62), length(t))';
%Ynon = mat2cell(non_obese_data, ones(1, 62), length(t))';
%BHM_non_obese = bhmcmc(Ynon, Tnon, 5, 2, 1/100, 1, 2000, 10000, 0, [], [], []); % for common grid

save('ObeseAug10.mat', 'BHM_obese', 'Sigma_est_ob', 'mu_p_ob')

% get ObeseJuly1.mat
% bhmcmc(Y, T, delta, cgrid, b, bs, Burnin, M, mat, Sigma_est, mu_p, nu)

       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       


  


