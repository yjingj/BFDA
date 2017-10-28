%% Author: Jingjing Yang (yjingj@gmail.com)
% Main function to implement the MCMC algorithm for smoothing and estimation of
%   functional data by BABF --- Bayesian process with approximation by cubic splines for random-grid data
%% Inputs
 % Y: 1 by n cell array with raw data, each element is an observation;
 % T: 1 by n cell array with grids for Y;
 % delta: determins the degrees of freedom of IW is delta+p-1; usually use 5
 % Burnin, M : number of iterations for MCMC
 % mat: whether use Matern kernel in IWP
 % Sigma_est, mu_est, tau: mean/covariance estimate on grid tau
 % c: determine the prior for the mean function
 % p: length of tau (work grid)
 % w, ws: decide hyper gamma priors
 % nu: order of smoothness in the Matern function; usually use 2.5
 % eval_grid: the common grid that mean-covariance functions will be
 % evaluated on.

%% Outputs a structure with elements
 % Zt: Bayesian Signals Estimates (smoothed) corresponding to Tcell, cell of 1 by n;
 % Zeta, Zeta_CL, Zeta_UL: Coefficients estimates, lower and upper 95% credible intervals;
 % Sigma_tau: Covariance matrix of functional data on working grid tau;
 % Sigma_zeta, Sigma_zeta_CL, Sigma_zeta_UL: Covariance matrix of basis function coefficients, lower and upper 95% credible intervals;
 % mu_tau: curve mean on working grid tau;
 % mu_zeta, nu_zeta_CI: mean of basis function coefficients, and 95% credible intervals;
 % Z_cgrid, Sigma_cgrid, mu_cgrid: signal, covariance, mean estiamtes on
 %                              the common eval_grid;
 % Z_cgrid_CL, Z_cgrid_UL, Sigma_cgrid_CL, Sigma_cgrid_UL, mu_cgrid_CI: 95% credible itervals;
 % rn, rn_CI: \gamma_{noise}, noise precision, and 95% credible iterval;
 % rs, rs_CI: \sigma_s^2, and 95% credible iterval;
 % rho, nu: if using matern covariance function
 % along with 95% credible intervals
 % Btau, BT: basis function evaluations on tau and T.
 % optknots: selected optimal knots by optknt() for constructing cubic B-splines
%% 

function [output] = babf_mcmc(Y, T, delta, Burnin, M, mat, Sigma_est, mu_est, tau, w, ws, c, nu, eval_grid)

%% working grid tau
P = cellfun(@length, Y); % length of all observation grids
n = size(Y, 2); % # of signals
m = length(tau); % # of working grids
J = ones(m, 1);
D = abs(J * tau - tau' * J'); % distance matrix on the working grid tau
p_cgrid = length(eval_grid);
display(['eval_grid has length ', p_cgrid])

%% MCMC sampling set up

%delta = 5; c = 1; %Choose prior normal distribution for mean (arbitrarily)
dm = delta + m - 1; %degrees of freedom for the inv-Wishart with delta
dmp = n + delta + m; %degrees of freedom for the conditional posterior inv-Wishart with n+1+delta

%% Solve for rho and nu
% by minimizing the MSE between A and Matern correlation matrix
if (mat) 
    
    A = COR(Sigma_est);
    options = optimoptions('fmincon','Display','notify');
    
    if(isempty(nu))       
        myobj = @(x) mean(mean((A - Matern(D, x(1), x(2), 1)).^2)); 
        [xx, fval] = fmincon(myobj, [1; 3], [], [], [], [], ...
            [1e-3; 2.5], [Inf; 4], [],options);
        rho = xx(1);
        nu = xx(2);  

        if(isnan(rho) )
        display('solve for rho and nu failed, fix nu = 2.5')
        nu = 2.5;
        [xx, fval] = fmincon(myobj, [1; nu], [], [], [], [], ...
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
  display([ '; rho = ', num2str(rho), '; nu = ', num2str(nu)]);

  A = Matern(D, rho, nu, 1); 
  
else
    A = Sigma_est;
    %A = topdm(Sigma_est);
    display('using empirical estimated covariance structure in the IW scale matrix');
    
end


%% calculate basis functions on tau and T, using bspline package
optknots = optknt(tau, 4); % select optimal knots sequences for working grid tau
Btau = bspline_basismatrix(4, optknots, tau);

BT_cgrid = bspline_basismatrix(4, optknots, eval_grid);

BT = cell(1, n);
for i = 1: n
    BT{i} = bspline_basismatrix(4, optknots, T{i});
end 


%% Empirical estimates of noise variance and \sigma_s^2

func = @(x) sum(diff(x) .^ 2);
snhat2 =  sum(cellfun(func, Y)) / (2 * sum(P - 1)); 
rn = 1 / snhat2; % Sample estimated noise precision
rs = (trace(Sigma_est) - m * snhat2 ) / (trace(A) / (delta - 2)); 

display(['Initial Estimates: ', ' noise variance = ', num2str(snhat2), ...
    '; sigma_s^2 = ', num2str(rs)])

%Choose inv-gamma distribution for noisy variance w = 1; ws = 5;
b = 1 / w; a = b * rn;
bs = 1 / ws; as = bs * rs;

%Choose gamma distribution for 1/noise-variance and signal-variance
display([' a = ',num2str(a),'; b = ',num2str(w),'; as = ',num2str(as),'; bs = ',num2str(ws)])


%% Setup initial values

mu_est = reshape(mu_est, m, 1);
mu0 = mu_est; 

rn = 1 / snhat2; % Sample estimated noise precision
rs = (trace(Sigma_est) - m * snhat2 ) / (trace(A) / (delta - 2)); 

Zt = cell(1, n); 
iK_tau = Sigma_est;

sumsquare_func = @(x, y) sum((x - y).^2);

% mean and covariance of coefficients zeta
mu_zeta = Btau \ mu_est; % p X 1 dimension
iK_zeta = (Btau \ iK_tau) / Btau'; % p by p dimension
K_zeta = pinv(iK_zeta); 

Zeta = repmat(mu_zeta, 1, n);
mu0_zeta = Btau \ mu0;
A_zeta = (Btau\A) / (Btau'); 

%% Assign memory ahead  to store MCMC samples

mu_cgrid_Out = NaN(p_cgrid, M); % Mean vector over eval_grid
iK_cgrid_Out = NaN(p_cgrid, p_cgrid, M); % Covariance matrix over eval_grid
Z_cgrid_Out = NaN(p_cgrid, n, M); % Signal mcmc samples over eval_grid

mu_zeta_Out = NaN(m, M); % Mean vector of zeta 
iK_zeta_Out = NaN(m, m, M); % Covariance matrix of zeta
ZetaOut = NaN(m, n, M); % Zeta mcmc samples
rnOut = NaN(1, M); % Precision of error, 1/noise-variance
rsOut = NaN(1, M); % signal variance, \sigma_s^2

%% Gibbs sampler

display('Start MCMC') 

for iter = 1 : (M + Burnin)    

    % update zeta/Z
        for i = 1:n
          
           iK_zeta_i = pinv(BT{i}' * BT{i} * rn + K_zeta);
           mu_zeta_i = iK_zeta_i * (BT{i}' * Y{i}' * rn + K_zeta * mu_zeta);
            
            Zeta(:, i) = mu_zeta_i + mychol(iK_zeta_i) * normrnd(0, 1, m, 1);
            Zt{i} = (BT{i} * Zeta(:, i))'; % approximate Z_i(ti)
        end
        
    %update 1/(noise variance), rn   
    rn = gamrnd(sum(P)/2 + a, 1/(b + sum(cellfun(sumsquare_func, Y, Zt)) / 2) );
    
    %update mu_zeta
    Mu_var = iK_zeta ./ (n + c);
    Mu_mean = (sum(Zeta, 2) + c .* mu0_zeta) ./ (n + c);    
    mu_zeta = Mu_mean + mychol(Mu_var) * normrnd(0, 1, m, 1);
    
    %update precision matrix K_zeta
    G = ((Zeta - repmat(mu_zeta, 1, n)) * (Zeta - repmat(mu_zeta, 1, n))' +...
        c .* (mu_zeta -  mu0_zeta) * (mu_zeta -  mu0_zeta)' + rs .* A_zeta);
    [iK_zeta, K_zeta] = myiwishrnd(dmp, G, m);
    
     %update signal variance rs
     iK_tau = (Btau) * iK_zeta * (Btau') ;
     K_tau = pinv(iK_tau) ;
     rs = gamrnd(m * dm/2 + as, 1/(bs + trace(A * K_tau) / 2));  

    % Save all MCMC samples   
    if iter > Burnin
       mu_zeta_Out(:, (iter - Burnin)) = mu_zeta;
       iK_zeta_Out(:, :, (iter - Burnin)) = iK_zeta; 
       ZetaOut(:, :, (iter - Burnin)) = Zeta;
       rnOut(iter-Burnin) = rn;
       rsOut(iter - Burnin) = rs;
       
       mu_cgrid_Out(:, (iter - Burnin)) = BT_cgrid * mu_zeta;
       iK_cgrid_Out(:, :, (iter - Burnin)) = BT_cgrid * iK_zeta * BT_cgrid';
       Z_cgrid_Out(:, :, (iter - Burnin)) = (BT_cgrid * Zeta);
       
    end
    
end

display('End of MCMC ... ')


%% MCMC diagnosis using mcmcdiag library
 display(['Calculate Potential Scale Reduction Factor (PSRF)...']);
 display('PSRF < 1.2 means the MCMC chain mixed well and acheived convergence.');
 
 display(['PSRF for 1/sigma_epsilon^2: ', num2str( psrf(rsOut'))]);
 display(['PSRF for 1/sigma_s^2: ', num2str( psrf(rnOut'))]);
 display(['PSRF for Zeta(1,1) : ', num2str(psrf(reshape(ZetaOut(1, 1, :), M, 1)))]);
 display(['PSRF for Sigma_zeta(1,1) : ', num2str(psrf(reshape(iK_zeta_Out(1, 1, :), M, 1)))]); 
 display(['PSRF for mu(1): ', num2str( psrf(mu_zeta_Out(1, :)'))]);
 
%% Calculate MCMC sample average
display('Calculating posterior sample means...');

mu_cgrid = mean(mu_cgrid_Out(:, :), 2);
iK_cgrid = mean(iK_cgrid_Out(:, :, :), 3);
Z_cgrid = mean(Z_cgrid_Out(:, :, :), 3);

Zeta = mean(ZetaOut(:, :, :), 3); 
iK_zeta = mean(iK_zeta_Out(:, :, :), 3);
mu_zeta = mean(mu_zeta_Out(:, :), 2);
rs = mean(rsOut);
rn = mean(rnOut);
iK_zeta_SE = cov(Zeta');

% Approximate functional signals, mean, covariance
for i = 1:n
    Zt{i} = (BT{i} * Zeta(:, i))'; 
end
mu_tau = Btau * mu_zeta;
iK_tau = Btau * iK_zeta * Btau';

%% 95% pointwise confidence interval
Zeta_sort = sort(ZetaOut, 3); 
iK_zeta_sort = sort(iK_zeta_Out, 3); 
mu_zeta_sort = sort(mu_zeta_Out, 2);

Z_cgrid_sort = sort(Z_cgrid_Out, 3); 
iK_cgrid_sort = sort(iK_cgrid_Out, 3); 
mu_cgrid_sort = sort(mu_cgrid_Out, 2);

rs_sort = sort(rsOut); 
rn_sort = sort(rnOut); 

q1 = ceil(0.025 * M);
q2 = floor(0.975 * M);

Zeta_CL = Zeta_sort(:, :, q1);
Zeta_UL = Zeta_sort(:, :, q2);
iK_zeta_CL = iK_zeta_sort(:, :, q1);
iK_zeta_UL = iK_zeta_sort(:, :, q2);
mu_zeta_CI = [mu_zeta_sort(:, q1), mu_zeta_sort(:, q2)];

Z_cgrid_CL = Z_cgrid_sort(:, :, q1);
Z_cgrid_UL = Z_cgrid_sort(:, :, q2);
iK_cgrid_CL = iK_cgrid_sort(:, :, q1);
iK_cgrid_UL = iK_cgrid_sort(:, :, q2);
mu_cgrid_CI = [mu_cgrid_sort(:, q1), mu_cgrid_sort(:, q2)];

rs_CI = [rs_sort(q1), rs_sort(q2)];
rn_CI = [rn_sort(q1), rn_sort(q2)];


%% outputs

if(mat)
    output = struct('Zt', {Zt}, 'Sigma_tau', iK_tau, 'Sigma_zeta_SE', iK_zeta_SE, ...
    'mu_tau', mu_tau, 'mu_zeta', mu_zeta, 'rn', rn, 'rs', rs,...
    'rho', rho, 'nu', nu, 'Zeta_CL', Zeta_CL, ...
    'Zeta_UL', Zeta_UL, 'Sigma_zeta_CL', iK_zeta_CL, ...
    'Sigma_zeta_UL', iK_zeta_UL, 'mu_zeta_CI', mu_zeta_CI, ...
    'Btau', Btau, 'BT', {BT}, 'Zeta', Zeta, 'Sigma_zeta', iK_zeta, ...
    'Z_cgrid', Z_cgrid, 'Z_cgrid_CL', Z_cgrid_CL, 'Z_cgrid_UL', Z_cgrid_UL, ...
    'Sigma_cgrid', iK_cgrid, 'Sigma_cgrid_CL', iK_cgrid_CL, 'Sigma_cgrid_UL', iK_cgrid_UL,...
    'mu_cgrid', mu_cgrid, 'mu_cgrid_CI', mu_cgrid_CI, ...
    'rs_CI', rs_CI, 'rn_CI', rn_CI, 'optknots', optknots);
else
    output = struct('Zt', {Zt}, 'Sigma_tau', iK_tau, 'Sigma_zeta_SE', iK_zeta_SE, ...
    'mu_tau', mu_tau, 'mu_zeta', mu_zeta, 'rn', rn, 'rs', rs,...
    'Zeta_CL', Zeta_CL, ...
    'Zeta_UL', Zeta_UL, 'Sigma_zeta_CL', iK_zeta_CL, ...
    'Sigma_zeta_UL', iK_zeta_UL, 'mu_zeta_CI', mu_zeta_CI, ...
    'Btau', Btau, 'BT', {BT}, 'Zeta', Zeta, 'Sigma_zeta', iK_zeta, ...
    'Z_cgrid', Z_cgrid, 'Z_cgrid_CL', Z_cgrid_CL, 'Z_cgrid_UL', Z_cgrid_UL, ...
    'Sigma_cgrid', iK_cgrid, 'Sigma_cgrid_CL', iK_cgrid_CL, 'Sigma_cgrid_UL', iK_cgrid_UL,...
    'mu_cgrid', mu_cgrid, 'mu_cgrid_CI', mu_cgrid_CI, ...
    'rs_CI', rs_CI, 'rn_CI', rn_CI, 'optknots', optknots);
end

display('BABF mcmc completed.');




