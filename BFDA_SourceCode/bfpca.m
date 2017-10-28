% Author: Jingjing Yang (yjingj@gmail.com)
% Main function to implement the MCMC algorithm for BFPCA
 
%% Inputs
 % Y: 1 by n cell array with raw data, each element is an observation;
 % T: 1 by n cell array with common grids;
 % Burnin, M : number of iterations for MCMC
 % Sigma_est: empirical covariance estimate
 % a, b, w: determin shape, scale parameters for IG of lambda
 
 %% Outputs a structure with elements
 % Z: smoothes signals

 % Zeta, Zeta_CL, Zeta_UL: Bayesian Estimates and 95% credible intervals for the PC coefficients, K by n;
 % rn, rn_CI, lambda, lambda_CI: Bayesian estimates and 95% credible intervals of the noise precision, PC coefficient precisions;
 % mu: sample estimate, p by 1;
 % Sigma_est: covariance estimate;

%% 
function [output] = bfpca(Y, T, Burnin, M, Sigma_est, a, b, w)


%%
t = T{1};
p = length(t); % length of the pooled grid
n = size(Y, 2); % # of signals

% J = ones(p, 1);
% D = abs(J * pgrid - pgrid' * J');

Ymat = reshape(cell2mat(Y), [p, n]); % noisy signal data


%% PCA of Sigma_est

[V, E] = eig(Sigma_est);
eig_val = diag(E);
K = find(cumsum(eig_val / sum(eig_val) ) > 0.99);
K = K(1); % select # of PCs

display(['Selected ', num2str(K), 'PCs to explain 99% variance.'])

%% assign memory ahead
% Setup empty data matrix to store MCMC samples
lambdaOut = NaN(K, M); % Mean vector
ZetaOut = NaN(K, n, M); % Smoothed curve
rnOut = NaN(1, M); % Precision of error, \gamma_n

%% Initial set up parameter values
Psi = V(:, 1:K); % eigen functions
mu = nanmean(Ymat, 2); % Sample mean
lambda = ones(K, 1); % precision for zeta

% Determine hyper-priors for noise precision
func = @(x) sum(diff(x) .^ 2);
snhat2 =  sum(cellfun(func, Y)) / (2 * n * (p - 1)); 
rn = 1 / snhat2; 
b_noise = 1 / w; a_noise = b_noise * rn;


   
%% Gibbs sampler

display('Starting MCMC...') 

for iter = 1 : (M + Burnin)    

    Vzeta = diag(1./(lambda + rn));
    Zeta_mean = (Vzeta) * (rn .* Psi' * Ymat ); 
    L = sqrt(Vzeta);
    Zeta = Zeta_mean + L * normrnd(0, 1, K, n); 

    %update 1/(noise variance), rn 
    Z = Psi * Zeta; % recover signals
    rn = gamrnd((n*p)/2 + a_noise, 1/(b_noise + nansum(nansum((Ymat - Z).^2))/2));

    %update lambda
    for k = 1:K
        lambda(k) = gamrnd(n/2 + a, 1/(b + sum(Zeta(k, :).^2) / 2)); 
    end
    
    % Save all MCMC samples   
    if iter > Burnin
       lambdaOut(:, (iter - Burnin)) = lambda;
       ZetaOut(:, :, (iter - Burnin)) = Zeta; 
       rnOut(iter - Burnin) = rn;
    end
    
end

display('Ending MCMC...')

%% MCMC diagnosis
 display(['Calculate Potential Scale Reduction Factor (PSRF)...']);
 display('PSRF < 1.2 means the MCMC chain mixed well and acheived convergence.');
 
 display(['PSRF for 1/sigma_s^2: ', num2str( psrf(rnOut'))]);
 display(['PSRF for Z(1,1) : ', num2str(psrf(reshape(ZetaOut(1, 1, :), M, 1)))]);
 display(['PSRF for lambda(1): ', num2str( psrf(lambdaOut(1, :)'))]);


%% Calculate MCMC sample average
display('Calculating posterior sample means...');
Zeta = mean(ZetaOut(:, :, :), 3); 
lambda = mean(lambdaOut(:, :), 2);
rn = mean(rnOut);
Z = Psi * Zeta;

%% 95% pointwise confidence interval
display('Calculating 95% CI...')
Zeta_sort = sort(ZetaOut, 3); 
lambda_sort = sort(lambdaOut, 2);
rn_sort = sort(rnOut); 

q1 = ceil(0.025 * M);
q2 = floor(0.975 * M);

Zeta_CL = Zeta_sort(:, :, q1);
Zeta_UL = Zeta_sort(:, :, q2);
lambda_CI = [lambda_sort(:, q1), lambda_sort(:, q2)];
rn_CI = [rn_sort(q1), rn_sort(q2)];


%%
output = struct('Z', Z, 'Zeta', Zeta, 'Zeta_CL', Zeta_CL, 'Zeta_UL', Zeta_UL, ...
        'rn', rn, 'rn_CI', rn_CI, ...
        'lambda', lambda, 'lambda_CI', lambda_CI, ...
        'mu', mu, 'Sigma_est', Sigma_est);
    
display('BFPCA completed.');

end