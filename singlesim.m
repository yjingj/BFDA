
 clear all, 
    
 %%
stream = RandStream('twister','Seed', 2013);
reset(stream);  % set up a seed for simulation

n = 100; % Number of curves
p = 30; % Numberof pooled grid points
t = (0:1/(p-1):1); %Pooled grid

snr = 2; % signal-to-noise ratio
ell = 1/2; sf = sqrt(5); sn = sf/ snr;
%ell is rho in our model
%sn is the true standard deviation of noise 1.1180
%sf is the true standard deviation of signal  2.2361

C = zeros(p); 
D = zeros(p); % Distance matrix for the pooled grid
for i = 1 : p
    for j = 1 : p
    D(i, j) = abs(t(i) - t(j));
    C(i, j) = Matern(D(i, j), ell, 2.5, sf^2); %degree of Matern Covariance Function is nu = 2.5  
    end
end
A0 = C ./ sf^2;  %True correlation matrix on pooled grid

%%
mut = 3*sin(4*t); % True signal mean on pooled grid
P = zeros(1, n); % # of grid points for each signal

Idx = cell(1, n); % Cell of idxes of observations on the pooled grid
T = cell(1, n); % Cell of grids
X = cell(1, n); % Cell of true data
Y = cell(1, n); % Cell of raw data
MU = cell(1, n); %Cell of means

%Generate functional data on pooled grid
X0 = chol(C)' * normrnd(0, 1, p, n) + repmat(mut', 1, n); % True signal, size p x n
Y0 = X0 + sn * normrnd(0, 1, p, n); % Raw signall = Truth + Noisy

%% Create functional data with the common pooled grid
for i = 1 : n
    Idx{i} = 1:p;%sort(randsample(1:p, 0.8 * p));
    T{i} = t(Idx{i}); 
    MU{i} = mut(Idx{i});
    P(i) = length(Idx{i});
    
    X{i} = X0(Idx{i}, i)'; %True signal
    Y{i} = Y0(Idx{i}, i)';    
end 

%% Create functional data with uncommon grid
for i = 1 : n
    Idx{i} = sort(randsample(1:p, 0.8 * p)); 
    % 80% grid points are randomly sampled from the pooled grid
    T{i} = t(Idx{i}); 
    MU{i} = mut(Idx{i});
    P(i) = length(Idx{i});
    
    X{i} = X0(Idx{i}, i)'; %True signal
    Y{i} = Y0(Idx{i}, i)';    
end 
%%
plot(Y0)


%% PACE
addpath(genpath('/Users/jingjingyang/Documents/MATLAB/BSCG/PACE_example/release2.11'))

param_X = setOptions('regular',0,'FVE_threshold',0.99,'corrPlot',0,'rho',-1); 
% regular=0 corresponds to the sparse case.   
 % regular=2 corresponds to the dense data.   
 
Simu_pcaRes= FPCA(Y, T, param_X);   %perform PCA on x

% Take the results out of the PCA analysis.
K_x = getVal(Simu_pcaRes,'no_opt'); % the number of Eigenfunctions chosen.
pcx = getVal(Simu_pcaRes,'xi_est');  % The estimated PC scores. 
Phihat=getVal(Simu_pcaRes,'phi');  % The estimated Eigenfunctions.
snhat2=getVal(Simu_pcaRes,'sigma'); % The estimated noise for X(t)
lamhat=getVal(Simu_pcaRes,'lambda'); % The estimated eigenvalues.

Est_Sigma=Phihat*diag(lamhat)*Phihat'; % The estimated covariance
R = Est_Sigma ./ (sqrt(diag(Est_Sigma)) * sqrt(diag(Est_Sigma))'); % Estiamted correlation

figure()
plot(Phihat)

[xgrid,ygrid]=meshgrid(t); % 3D plots of the covariance matrices
figure()
mesh(xgrid,ygrid,R)
rmse(A0, R)

%% Bayesian Smoothing by MCMC algorithm 

[Z, iK, mu, rn, rs, rt] = bsmoothos(t, Y0, 2.5, 5); % for common grid
% bsmoothos function is to be used for common grid


[Z, iK, mu, rn, rs, rt, Est_Sigma, t] = bhmcmc(Y, T, 2.5, 5, 0);
save('sparsesim.mat', 'Z', 'iK', 'mu', 'rn', 'rs', 'rt', 'Est_Sigma', 't', 'Y', 'Y0', 'X0', 'A0', 'C', 'T', 'Idx', 'Idx_rest', 'Yfull')


% delta = 5;
% nu=2.5;
% cgrid = 0 denotes uncommon grid

