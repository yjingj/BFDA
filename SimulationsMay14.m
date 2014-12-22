%clear all;
stream = RandStream('twister','Seed', 2014);
reset(stream);  % set up a seed for simulation

%% Parameter set up
addpath(genpath(cat(2, pwd, '/PACErelease2.11')))
n = 100; % Number of curves
p = 40; % Numberof pooled grid points
t = (0 : (pi/2)/(p-1) : (pi/2)); %Pooled grid

[xgrid,ygrid]=meshgrid(t); % 3D plots of the covariance matrices

sf = sqrt(5); % singal standard deviation
snr = 2; % signal to noise ratio
rho = 1/2; % spacial scale parameter
nu = 2.5;
dense = 0.8; % percentage of data density on pooled grid
% Repeat the simulation for 100 times

Burnin = 2000; M = 10000;
b = 1; bs = 5;  delta = 30;  cgrid = 2;

param_X = setOptions('regular', 2,'FVE_threshold',0.99,'corrPlot',0,'rho',-1); 
param_X_sparse = setOptions('regular', 0, 'FVE_threshold', 0.99, 'corrPlot', 0, 'rho', -1); 


%%
Nsim = 100;
%YSim = zeros(p, n, Nsim);
%XSim = zeros(p, n, Nsim);

YSim_NS = zeros(p, n, Nsim);
XSim_NS = zeros(p, n, Nsim);

% save results form BHM
%ZSim = zeros(p, n, Nsim); % Smoothed curve
%iKSim = zeros(p, p, Nsim); % Covariance matrix
%muSim = zeros(p, Nsim); % Mean vector
%rnSim = zeros(1, Nsim); % Precision of error, gamma_n
%rsSim = zeros(1, Nsim); % signal variance, sf
%rhoSim = zeros(1, Nsim); % gamma_theta1

% save results from BHM_sparse
ZSim_NS = zeros(p, n, Nsim); % Smoothed curve
iKSim_NS = zeros(p, p, Nsim); % Covariance matrix
muSim_NS = zeros(p, Nsim); % Mean vector
rnSim_NS = zeros(1, Nsim); % Precision of error, gamma_n
rsSim_NS = zeros(1, Nsim); % signal variance, sf
rhoSim_NS = zeros(1, Nsim); % gamma_theta1

% save results from PCA
%ZpSim = zeros(p, n, Nsim); % Smoothed curve
%iKpSim = zeros(p, p, Nsim); % Covariance matrix
%mupSim = zeros(p, Nsim); % Mean vector
%snpSim = zeros(1, Nsim); % noise variance

ZpSim_NS = zeros(p, n, Nsim); % Smoothed curve
iKpSim_NS = zeros(p, p, Nsim); % Covariance matrix
mupSim_NS = zeros(p, Nsim); % Mean vector
snpSim_NS = zeros(1, Nsim); % noise variance


%%  

for (iSim = 1 : Nsim)

    %************************** Generate stationary Gaussian data
    %GausData = gaussian_sim(t, n, sf, snr, nu, rho, dense);
    %Y = GausData.Y_cgrid;
    %T = GausData.T_cgrid;
    
    % simulated true and noise functional signals for non-sparse case
   % XSim(:, :, iSim) = GausData.X0;
    %YSim(:, :, iSim) = GausData.Y0;
    
    %BHM_cgrid = bhmcmc_fix_rho(Y, T, nu, delta, cgrid, b, bs, Burnin, M);  % BHM model  
    %ZSim(:, :, iSim) = BHM_cgrid.Z; % Smoothed curve
    %iKSim(:, :, iSim) = BHM_cgrid.iK; % Covariance matrix
    %muSim(:, iSim) = BHM_cgrid.mu; % Mean vector 
    %rnSim(:, iSim) = BHM_cgrid.rn; % Precision of error, gamma_n
    %rsSim(:, iSim) = BHM_cgrid.rs; % signal variance, sf
    %rhoSim(:, iSim) = BHM_cgrid.rho; % gamma_theta1
    
    % PACE for non-stationary case
    %PCA_cgrid= FPCA(Y, T, param_X);   %perform PCA on x, common grid case
    %Phihat=getVal(PCA_cgrid,'phi');  % The estimated Eigenfunctions.
    %lamhat=getVal(PCA_cgrid,'lambda'); % The estimated eigenvalues.
    %pcx = getVal(PCA_cgrid,'xi_est');  % The estimated PC scores. 

    %mupSim(:, iSim) = (getVal(PCA_cgrid, 'mu'))';
    %snpSim(iSim) = getVal(PCA_cgrid,'sigma'); 
    %ZpSim(:, :, iSim) = (pcx * Phihat')';
    %iKpSim(:, :, iSim) = Phihat*diag(lamhat)*Phihat'; % The estimated covariance

     %************************ Genrate Non-Stationary Gaussian data
    GausData_NS = gaussian_simNS(t, n, sf, snr, nu, rho, dense);
    Y_NS = GausData_NS.Y_cgrid;
    T_NS = GausData_NS.T_cgrid;
    
    XSim_NS(:, :, iSim) = GausData_NS.X0;
    YSim_NS(:, :, iSim) = GausData_NS.Y0;
    
    %BHM_cgrid
    BHM_NS = bhmcmc_fix_rho(Y_NS, T_NS, nu, delta, cgrid, b, bs, Burnin, M);  % BHM model
    
    ZSim_NS(:, :, iSim) = BHM_NS.Z; % Smoothed curve
    iKSim_NS(:, :, iSim) = BHM_NS.iK; % Covariance matrix
    muSim_NS(:, iSim) = BHM_NS.mu; % Mean vector 
    rnSim_NS(:, iSim) = BHM_NS.rn; % Precision of error, gamma_n
    rsSim_NS(:, iSim) = BHM_NS.rs; % signal variance, sf
    rhoSim_NS(:, iSim) = BHM_NS.rho; % gamma_theta1
    
    % PACE 
    PCA_NS= FPCA(Y_NS, T_NS, param_X);   %perform PCA on x, common grid case
    Phihat=getVal(PCA_NS,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_NS,'lambda'); % The estimated eigenvalues.
    pcx = getVal(PCA_NS,'xi_est');  % The estimated PC scores. 

    mupSim_NS(:, iSim) = (getVal(PCA_NS, 'mu'))';
    snpSim_NS(iSim) = getVal(PCA_NS,'sigma'); 
    ZpSim_NS(:, :, iSim) = (pcx * Phihat')';
    iKpSim_NS(:, :, iSim) = Phihat*diag(lamhat)*Phihat'; % The estimated covariance

end

C = GausData.C;
mut = GausData.mut;
D = GausData.D;

C_NS = GausData_NS.C;
mut_NS = GausData_NS.mut;

save('StatSimMay14.mat', 'XSim', 'YSim', 'XSim_NS', 'YSim_NS', ...
     'ZSim', 'iKSim', 'muSim', 'rnSim', 'rsSim', 'rhoSim', ...
    'ZSim_NS', 'iKSim_NS', 'muSim_NS', 'rnSim_NS', 'rsSim_NS', 'rhoSim_NS', ...
    'ZpSim', 'snpSim', 'iKpSim', 'mupSim', 'ZpSim_NS', 'snpSim_NS', 'iKpSim_NS', 'mupSim_NS',...
    'C', 'mut', 'D', 'C_NS', 'mut_NS');

%% load('Device1_WRL.mat') real spectroscopy data from el paso 2013
%wcal 291 by 462
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

nu=2.5; delta = 5; cgrid = 2; b = 1; bs = 5; Burnin = 2000; M = 10000;
spec_BHM = bhmcmc_fix_rho(Y, T, nu, delta, cgrid, b, bs, Burnin, M, 2);

%%
% PACE 
addpath(genpath(cat(2, pwd, '/PACErelease2.11'))) 
param_X = setOptions('regular', cgrid,'FVE_threshold',1,'corrPlot',0,'rho',-1); 

    PCA_spec= FPCA(Y, T, param_X);   %perform PCA on x, common grid case
    Phihat=getVal(PCA_spec,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_spec,'lambda'); % The estimated eigenvalues.
    pcx = getVal(PCA_spec,'xi_est');  % The estimated PC scores. 

    mu_spec_pace = (getVal(PCA_spec, 'mu'))';
    sn_spec_pace = getVal(PCA_spec,'sigma'); 
    Z_spec_pace = (pcx * Phihat')' + repmat(mu_spec_pace, 1, 462);
    iK_spec_pace = Phihat*diag(lamhat)*Phihat'; % The estimated covariance

 save('emWL', 'Y', 'T', 'x', 'Yspecmat', 'spec_BHM', 'Z_spec_pace', 'mu_spec_pace', ...
           'iK_spec_pace', 'sn_spec_pace',  'PCA_spec',...
           'emWL', 'Yspecleft1', 'xleft1', 'Yspecleft2', 'xleft2')
       

       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       


  


