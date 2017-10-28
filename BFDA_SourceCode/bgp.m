% Author: Jingjing Yang (yjingj@gmail.com)
% Main function to implement the MCMC algorithm for Bayesian GP model
% Only for common grid case
 
%% Inputs
 % Y: 1 by n cell array with raw data, each element is an observation;
 % T: 1 by n cell array with grids for Y;
 % Burnin, M : number of iterations for MCMC
 % mat: stationary or not
 % Sigma_est, mu_est: empirical covariance/mean estimate
 % c: scale of prior mean variance 
 % w, ws: scale the variance for the prior of sigma_s^2
 % nu: given order for matern function
 % tol: singular values < tol were set to 0 in pinv()
  
 %% Outputs a structure with elements
 % Z, Z_CL, Z_UL: Bayesian Signals Estimates and 95% credible interval (smoothed), p by n;
 % mu, mu_CI: Bayesian estimate of curve mean, p by 1;
 % rn, rn_CI: Bayesian estimate and 95% credible interval of noise precision \gamma_{noise};
 % lambda, lambda_CI: Bayesian estiamte and 95% credible interval
 % rho, nu: Empirical estimate of the spacial scale parameter and order if
 % Sigma_est: Input covariance estimate, p by p;

%%
function [output] = bgp(Y, T, Burnin, M, mat, Sigma_est, mu_est, c, w, ws, nu, tol)
               
pgrid = sort(unique(cell2mat(T))); %pooled grid
p = length(pgrid); % length of the pooled grid
I = eye(p); % identity matrix
n = size(Y, 2); % # of signals

Yfull = reshape(cell2mat(Y), p, n);
J = ones(p, 1);
D = abs(J * pgrid - pgrid' * J');
     
%% Solve for rho and nu
% by minimizing the MSE between A and Matern correlation matrix

if (mat) 
    A = COR(Sigma_est);
    options = optimoptions('fmincon','Display','notify');
    
    if(isempty(nu))       
        myobj = @(x) mean(mean((A - Matern(D, x(1), x(2), 1)).^2)); 
        [xx, fval] = fmincon(myobj, [1; 2.5], [], [], [], [], ...
            [1e-3; 2.5], [Inf; 10], [],options);
        rho = xx(1);
        nu = xx(2);  

        if(isnan(rho) )
        display('solve for rho and nu failed, fix nu = 2.5')
        nu = 2.5;
        [xx, fval] = fmincon(myobj, [1; nu], [], [], [], [],...
            [1e-3; nu], [Inf; nu], [],options);
        rho = xx(1);
        end

    else 
         myobj = @(x) mean(mean((A - Matern(D, x(1), x(2), 1)).^2)); 
         display('solve for rho with given value for nu')
        [xx, fval] = fmincon(myobj, [1; nu], [], [], [], [], ...
            [1e-3; nu], [Inf; nu], [],options);
        rho = xx(1);
    end
    
  display('using Matern correlation structure in the IW scale matrix');   
  display(['rho = ', num2str(rho), '; nu = ', num2str(nu)]);
  A = Matern(D, rho, nu, 1); 

else
    A = Sigma_est; %A = topdm(Sigma_est);
    display('using empirical estimated covariance structure in the IW scale matrix');
    
end


%%
% Empirical estimates of noise variance and \sigma_s^2
func = @(x) sum(diff(x) .^ 2);
snhat2 =  sum(cellfun(func, Y)) / (2 * n * (p - 1)); 
rn = 1 / snhat2;

lambda = trace(A) / (trace(Sigma_est) - p * snhat2 );

display(['Initial Estimates: ', ' noise variance = ', num2str(snhat2),'; sigma_s^2 = ', num2str(1/lambda)]);

% Determine hyper-priors
b = 1 / w; a = b * rn;
bs = 1 / ws; as = bs * lambda;

display(['Hyper-prior parameters: ', ' a = ', num2str(a), '; b = ', num2str(b),...
    '; as = ', num2str(as), '; bs = ', num2str(bs)])

% Initial set up parameter values
mu = nanmean(Yfull, 2); % Sample mean
mu0 = reshape(mu_est, p, 1);

iA = pinv(A, tol);
Sigma = A ./ lambda;
K = iA .* lambda; %pinv(Sigma);


%% assign memory ahead
% Setup empty data matrix to store MCMC samples
muOut = NaN(p, M); % Mean vector
ZOut = NaN(p, n, M); % Smoothed curve
rnOut = NaN(1, M); % Precision of error, \gamma_n
lambdaOut = NaN(1, M); % \sigma_s^2

   
%% Gibbs sampler

display('Starting MCMC...') 

for iter = 1 : (M + Burnin)    

         iZvar = (K + rn .* I);
         Zvar = pinv(iZvar, tol);
         Zmean = (iZvar) \ (rn .* Yfull + repmat(K * mu, 1, n)); 
         L = mychol(Zvar);
         Z = Zmean + L * normrnd(0, 1, p, n); 
         
    %update 1/(noise variance), rn   
    rn = gamrnd((n*p)/2 + a, 1/(b + nansum(nansum((Yfull - Z).^2))/2));
    
    %update signal mean mu
    Mu_var = Sigma ./ (n + c);
    Mu_mean = (sum(Z, 2) + c .* mu0) ./ (n + c); 
    L = mychol(Mu_var);
    mu = Mu_mean + L * normrnd(0, 1, p, 1);
   
    %update signal precision lambda 
    %lambda = gamrnd((p*n)/2 + as + p, ...
    %    1/(bs + trace(((Z- repmat(mu, 1, n))' /A) * (Z - repmat(mu, 1, n))) / 2 +...
    %    trace(c .* ((mu - mu0)' / A) * (mu - mu0))/2)); 


    lambda = gamrnd((p*n)/2 + as + p/2, ...
        1/(bs + trace((Z- repmat(mu, 1, n))' * iA * (Z - repmat(mu, 1, n))) / 2 +...
        c * trace((mu - mu0)' * iA * (mu - mu0))/2));
    
    %update signal precision matrix K
    Sigma = A ./ lambda;
    K = iA .* lambda; % pinv(Sigma); 
    
    % Save all MCMC samples   
    if iter > Burnin
       muOut(:, (iter - Burnin)) = mu;
       ZOut(:, :, (iter - Burnin)) = Z;
       rnOut(iter - Burnin) = rn;
       lambdaOut(iter - Burnin) = lambda;
    end
    
end

display('Ending MCMC...')



%% MCMC diagnosis
 display(['Calculate Potential Scale Reduction Factor (PSRF)...']);
 display('PSRF < 1.2 means the MCMC chain mixed well and acheived convergence.');
 
 display(['PSRF for lambda: ', num2str( psrf(lambdaOut'))]);
 display(['PSRF for noise precision: ', num2str( psrf(rnOut'))]);
 display(['PSRF for Z(1,1) : ', num2str(psrf(reshape(ZOut(1, 1, :), M, 1)))]);
 display(['PSRF for mu(1): ', num2str( psrf(muOut(1, :)'))]);


%% Calculate MCMC sample average
display('Calculating posterior sample means...');
Z = mean(ZOut(:, :, :), 3); 
mu = mean(muOut(:, :), 2);
lambda = mean(lambdaOut);
rn = mean(rnOut);

%% 95% pointwise confidence interval
display('Calculating 95% CI...')

Z_sort = sort(ZOut, 3); 
mu_sort = sort(muOut, 2);
lambda_sort = sort(lambdaOut); 
rn_sort = sort(rnOut); 


q1 = ceil(0.025 * M);
q2 = floor(0.975 * M);

Z_CL = Z_sort(:, :, q1);
Z_UL = Z_sort(:, :, q2);
mu_CI = [mu_sort(:, q1), mu_sort(:, q2)];
lambda_CI = [lambda_sort(q1), lambda_sort(q2)];
rn_CI = [rn_sort(q1), rn_sort(q2)];


%%
if(mat)
output = struct('Z', Z, 'Z_CL', Z_CL, 'Z_UL', Z_UL, ...
        'mu', mu, 'mu_CI', mu_CI, ...
        'rn', rn, 'rn_CI', rn_CI,...
        'lambda', lambda, 'lambda_CI', lambda_CI, ...
        'rho', rho, 'nu', nu);
else
    output = struct('Z', Z, 'Z_CL', Z_CL, 'Z_UL', Z_UL, ...
        'mu', mu, 'mu_CI', mu_CI, ...
        'rn', rn, 'rn_CI', rn_CI,...
        'lambda', lambda, 'lambda_CI', lambda_CI, ...
        'Sigma_est', Sigma_est);
end


display('BGP completed.');

end