function [output] = bhmcmc(Y, T, delta, cgrid, Burnin, M, mat, Sigma_est, mu_p, nu, ws, w)

% Author: Jingjing Yang (yjingj@gmail.com)
% Main function to implement the MCMC algorithm for smoothing functional
% data
 
%% Inputs
 % Y: 1 by n cell array with raw data, each element is an observation;
 % T: 1 by n cell array with grids for Y;
 % delta: determins the degrees of freedom of IW is delta+p-1; usually use 5
 % cgrid = 2 denote common grid case; while cgrid =0 denotes uncommon grid
 % b: decide hyper gamma priors
 % Burnin, M : number of iterations for MCMC
 % mat: 2 menas matern sturcture for Psi, 0 uses PACE correlation estimate
 % Sigma_est: empirical covariance estimate
 % mu_p: empirical mean estimate
 % A: empirical correlation estimate
 % nu: order of smoothness in the Matern function; usually use 2.5
 % ws: scale the variance for the prior of sigma_s^2
 % w: scale the variance for the prior of rn
 
 %% Outputs a structure with elements
 % t: pooled grid
 % Z: Bayesian Signals Estimates (smoothed), p by n;
 % iK: Bayesian estimate of Covariance matrix, p by p;
 % mu: Bayesian estimate of curve mean, p by 1;
 % rn: Bayesian estimate of noise precision \gamma_{noise};
 % rs: Bayesian estimate of the scala in the scale matrix of IW distribution \sigma_s^2;
 % rho, nu: Empirical estimate of the spacial scale parameter and order if
 % mat is true
 % parameter in Matern function, if not given 
 % Sigma_est: empirical covariance estimate by PACE, if not given
 % Yfull: raw data in matrix form
 % mat: given by input
 % along with 95% credible intervals for all Bayesian estimates
 
 %%                   
t = sort(unique(cell2mat(T))); %pooled grid
p = length(t); % length of the pooled grid
I = eye(p); % identity matrix
n = size(Y, 2); % # of signals
P = cellfun(@length, Y); % # of observations for each signal

D = zeros(p); % Distance matrix
for i = 1 : p
    for j = 1 : p
    D(i, j) = abs(t(i) - t(j));
    end
end

%% construct matrix form of the raw signal data

if  ~cgrid
    display('sparse case with uncommon grid')
    Idx = cell(n, 1); % Cell of index for observed grid on the pooled grid
    Idx_rest = cell(n, 1); % Cell of index for unobserved grid on the pooled grid
    Yfull = NaN(p, n); % n by p data matrix with nan's for unobserved data
    for i = 1:n
       Idx{i} = find(ismember(t, T{i}));
       Idx_rest{i} = find(~ismember(t, T{i}));
       Yfull(Idx{i}, i) = Y{i};
    end 
    
else 
    display('functional data with common grid')
    Yfull = reshape(cell2mat(Y), p, n);
end

%% MCMC sampling set up
%delta = 5
m = delta + p - 1; %degrees of freedom for the inv-Wishart with delta
mp = n + delta + p; %degrees of freedom for the conditional posterior inv-Wishart with n+1+delta
c = 1; %Choose prior normal distribution for mean (arbitrarily)

%% Find an empirical estimate for the covariance and sample mean by PACE, if not given

% Part of the 2D smoothing code from PACE. 
if (isempty(Sigma_est) || isempty (mu_p))
    addpath(genpath(cat(2, pwd, '/PACErelease2.11')))
    tt = cell2mat(T);  % 1 x N vector to hold the observed time points from all subjects
    yy = cell2mat(Y);  % 1 x N vector to hold the observed measurements from all subjects

    %Initial out1 is based on the unique time points of the pooled data.
    out1 =unique(tt);ngrid1 = length(out1); regular = cgrid; 
    error = 1; kernel = 'gauss'; verbose = 'on';

    % Determine bandwidth of mu using local smoothing using GCV.
    bw_mu = gcv_lwls(yy,tt,kernel,1,1,0,regular,verbose);    %use GCV method to choose bw for mean function
    bw_mu = adjustBW1(kernel,bw_mu,1,0,regular,verbose);

    %define the vector of case weight in the local weighted least square
    %here, it is set to be one for all subjects
    win1 = ones(1,length(tt));
    [~, mu_p] = lwls(bw_mu, kernel, 1,1,0,tt,yy',win1,out1); % mu is the locally smoothed mu(t) function.  
    rcov = getRawCov(Y, T, out1, mu_p, regular, 0); % Now get the raw Covariance, using (Y_{ij}-mu(t_{ij}))*(Y_{ij'}-mu(t_{ij'}))

    % Now smooth the raw covariance to get R estimate.
    [bw_xcov] =  gcv_mullwlsn(T,ngrid1,regular,error, kernel,rcov,verbose);
    bw_xcov = adjustBW2(kernel,bw_xcov,1,0,regular,verbose);

    [~,xcov]=mullwlsk(bw_xcov,kernel,rcov.tpairn,rcov.cxxn',rcov.win,t,t);  %smooth raw covariance;
    Sigma_est = (xcov+xcov')/2;   %transform the smoothed covariance matrix to guarantee it is a symmetric matrix.
    
else display('Run with given mean/covariance estimate ... ')

end


     
%% Solve for rho and nu
% by minimizing the MSE between A and Matern correlation matrix

if (mat) 
    
    A = COR(Sigma_est);
    
    if(isempty(nu))       
        myobj = @(x) mean(mean((A - Matern(D, x(1), x(2), 1)).^2)); 
        [xx, fval] = fmincon(myobj, [1; 2.5], [], [], [], [], [1e-3; 2.5], [Inf; 4.5]);
        rho = xx(1);
        nu = xx(2);  

        if(isnan(rho) )
        display('solve for rho and nu failed, fix nu = 2.5')
        nu = 2.5;
        [xx, fval] = fmincon(myobj, [1; nu], [], [], [], [], [1e-3; nu], [Inf; nu]);
        rho = xx(1);
        end

    else 
         myobj = @(x) mean(mean((A - Matern(D, x(1), x(2), 1)).^2)); 
         display('solve for rho with given value for nu')
        [xx, fval] = fmincon(myobj, [1; nu], [], [], [], [], [1e-3; nu], [Inf; nu]);
        rho = xx(1);
    end
    
  display('using Matern correlation structure in the IW scale matrix');  
  A = Matern(D, rho, nu, 1); 
  display(['minimum L2 distance', num2str(fval), '; rho = ', num2str(rho), '; nu = ', num2str(nu)]);

else
    
    A = Sigma_est;
    display('using empirical estimated covariance structure in the IW scale matrix');
    
end

%%
% Empirical estimates of noise variance and \sigma_s^2
func = @(x) sum(diff(x) .^ 2);
snhat2 =  sum(cellfun(func, Y)) / (2 * sum(P - 1)); 
rn = 1 / snhat2; 
rs = (trace(Sigma_est) - p * snhat2 ) / (trace(A) / (delta - 2)); 
display(['Initial Estimates: ', ' noise variance = ', num2str(snhat2),'; sigma_s^2 = ', num2str(rs)])

% Determine hyper-priors
b = rn / w;
a = b * rn;
%rs = ws * rs;
%display(['Hyper-prior parameters: ', ' a = ', num2str(a), '; b = ', num2str(b),'; rs = ', num2str(rs)])

bs = rn / ws;
as = bs * rs;
%bs = rs * ws;
%as = bs * rs * ws; 

display(['Hyper-prior parameters: ', ' a = ', num2str(a), '; b = ', num2str(b),...
    '; as = ', num2str(as), '; bs = ', num2str(bs)])

% Initial set up parameter values
mu = nanmean(Yfull, 2); % Sample mean
Z = Yfull;
Mu = repmat(mu, 1, n);
Z(isnan(Z)) = Mu(isnan(Z));

if(isempty(mu_p))
    mu0 = smooth(mu, 0.5, 'lowess');
else
    mu0 = mu_p';
end

iK = I;   % covariance matrix
K = pinv(iK); % precision matrix

%% assign memory ahead

% Setup empty data matrix to store MCMC samples
muOut = NaN(p, M); % Mean vector
iKOut = NaN(p, p, M); % Covariance matrix
ZOut = NaN(p, n, M); % Smoothed curve
rnOut = NaN(1, M); % Precision of error, \gamma_n
rsOut = NaN(1, M); % \sigma_s^2
   
%% Gibbs sampler

display('Starting MCMC...') 

for iter = 1 : (M + Burnin)    

    if ~cgrid       
        for i = 1:n          
            %mu3, V3 conditional mean/variance of Zi* | Zi, mu, \Sigma         
            %mu2, V2 conditional mean/variance of Zi | Zi*, mu, Yi \Sigma
            V11 = iK(Idx{i}, Idx{i});
            V12 = iK(Idx{i}, Idx_rest{i});
            V22 = iK(Idx_rest{i}, Idx_rest{i}); 
            
            B0 = V12' / V11;
            B1 = B0 * mu(Idx{i}) - mu(Idx_rest{i});
            mu3 = B0 * Z(Idx{i}, i) - B1;
            V3 = V22 - V12' / (V11) * V12;
            
            s=svd(V3);
            if s(1) < 0.0001
                Z(Idx_rest{i}, i) = mu3;
                iV2 = rn * eye(P(i)) + pinv(V11);
                V2 = pinv(iV2); 
                mu2 = iV2 \ (rn * Y{i}' + (V11) \ mu(Idx{i}));
            else
                L = mychol(V3);
                Z(Idx_rest{i}, i) = mu3 + L * normrnd(0, 1, length(Idx_rest{i}), 1);                
                B2 = B0'/ V3;
                iV2 = rn * eye(P(i)) + pinv(V11) + B2 * B0;
                V2 = pinv(iV2);
                mu2 = iV2 \ (rn * Y{i}' + V11 \ mu(Idx{i}) + B2 * (B1 + Z(Idx_rest{i}, i)));
            end
            
            L = mychol(V2);
            Z(Idx{i}, i) = mu2 + L * normrnd(0, 1, length(Idx{i}), 1);
        end
    else    
         iZvar = (K + rn .* I);
         Zvar = pinv(iZvar);
         Zmean = (iZvar) \ (rn .* Yfull + repmat(K * mu, 1, n)); 
         L = mychol(Zvar);
         Z = Zmean + L * normrnd(0, 1, p, n); 
    end

    %update 1/(noise variance), rn   
    rn = gamrnd(sum(P)/2 + a, 1/(b + nansum(nansum((Yfull - Z).^2))/2));
    
    %update signal mean mu
    Mu_var = iK ./ (n + c);
    Mu_mean = (sum(Z, 2) + c .* mu0) ./ (n + c); 
    L = mychol(Mu_var);
    mu = Mu_mean + L * normrnd(0, 1, p, 1);
    
    %update signal precision matrix K
    % Phi = rs .* A; 
    G = (Z- repmat(mu, 1, n)) * (Z - repmat(mu, 1, n))' +...
        c .* (mu - mu0) * (mu - mu0)' + rs .* A; 
    [iK, K] = myiwishrnd(mp, G, p);       
    
    %update signal variance rs 
    rs = gamrnd(p*m/2 + as, 1/(bs + trace(A * K) / 2)); 
    
    % Save all MCMC samples   
    if iter > Burnin
       muOut(:, (iter - Burnin)) = mu;
       iKOut(:, :, (iter - Burnin)) = iK; 
       ZOut(:, :, (iter - Burnin)) = Z;
       rnOut(iter - Burnin) = rn;
       rsOut(iter - Burnin) = rs;
    end
    
end

display('Ending MCMC...')

%% MCMC diagnosis
 addpath(genpath(cat(2, pwd, '/mcmcdiag')))
 display(['Calculate Potential Scale Reduction Factor (PSRF)...']);
 display('PSRF < 1.2 means the MCMC chain mixed well and acheived convergence.');
 
 display(['PSRF for 1/sigma_epsilon^2: ', num2str( psrf(rsOut'))]);
 display(['PSRF for 1/sigma_s^2: ', num2str( psrf(rnOut'))]);
 display(['PSRF for Z(1,1) : ', num2str(psrf(reshape(ZOut(1, 1, :), M, 1)))]);
 display(['PSRF for Sigma(1,1) : ', num2str(psrf(reshape(iKOut(1, 1, :), M, 1)))]); 
 display(['PSRF for mu(1): ', num2str( psrf(muOut(1, :)'))]);


%% Calculate MCMC sample average
display('Calculating posterior sample means...');
Z = mean(ZOut(:, :, :), 3); 
iK = mean(iKOut(:, :, :), 3);
mu = mean(muOut(:, :), 2);
rs = mean(rsOut);
rn = mean(rnOut);
iKSE = cov(Z');

%% 95% pointwise confidence interval
display('Calculating 95% CI...')
Z_sort = sort(ZOut, 3); 
iK_sort = sort(iKOut, 3); 
mu_sort = sort(muOut, 2);

rs_sort = sort(rsOut); 
rn_sort = sort(rnOut); 


q1 = ceil(0.025 * M);
q2 = floor(0.975 * M);

Z_CL = Z_sort(:, :, q1);
Z_UL = Z_sort(:, :, q2);
iK_CL = iK_sort(:, :, q1);
iK_UL = iK_sort(:, :, q2);
mu_CI = [mu_sort(:, q1), mu_sort(:, q2)];

rs_CI = [rs_sort(q1), rs_sort(q2)];
rn_CI = [rn_sort(q1), rn_sort(q2)];


%%
if(mat)
    output = struct('t', t, 'Z', Z, 'iK', iK, 'iKSE', iKSE, 'mu', mu, 'mu0', mu0, 'rn', rn, 'rs', rs,...
        'rho', rho, 'nu', nu, 'Sigma_est', Sigma_est, 'Z_CL', Z_CL, ...
        'Z_UL', Z_UL, 'iK_CL', iK_CL, 'iK_UL', iK_UL, 'mu_CI', mu_CI, ...
        'rs_CI', rs_CI, 'rn_CI', rn_CI,'Yfull', Yfull, 'mat', mat);
else
    output = struct('t', t, 'Z', Z, 'iK', iK, 'iKSE', iKSE, 'mu', mu, 'mu0', mu0, 'rn', rn, 'rs', rs,...
        'Sigma_est', Sigma_est, 'Z_CL', Z_CL, ...
        'Z_UL', Z_UL, 'iK_CL', iK_CL, 'iK_UL', iK_UL, 'mu_CI', mu_CI, ...
        'rs_CI', rs_CI, 'rn_CI', rn_CI,'Yfull', Yfull, 'mat', mat);
end
display('Job completed.');

end