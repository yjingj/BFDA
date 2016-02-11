function [output] = babf_mcmc(Y, T, delta, Burnin, M, mat, Sigma_tau, mu_tau_est, tau, p, w, ws, c, nu, reg, t_cgrid)

% To incooperate functional data with uncommon grid, sparse data
 
%% Inputs
 % Y: 1 by n cell array with raw data, each element is an observation;
 % T: 1 by n cell array with grids for Y;
 % delta: determins the degrees of freedom of IW is delta+p-1; usually use 5
 % Burnin, M : number of iterations for MCMC
 % mat: whether use Matern kernel in IWP
 % Sigma_tau, mu_tau, tau: mean/covariance estimate on grid tau
 % p: length of tau (work grid)
 % w, ws: decide hyper gamma priors
 % reg: regular code for pace
 % nu: order of smoothness in the Matern function; usually use 2.5
 % t_cgrid: the common grid that mean-covariance functions will be
 % evaluated on.
 
 
 %% Outputs a structure with elements
 % Zt: Bayesian Signals Estimates (smoothed), cell of 1 by n;
 % Zeta: Coefficients estimates, p by n;
 % iK_tau, iK_zeta: Covariance matrix, p by p;
 % mu_tau, mu_zeta: curve mean on working grid tau, p by 1;
 % rn: \gamma_{noise}, noise precision;
 % rs: \sigma_s^2;
 % rho, nu: if using matern covariance function
 % Est_Sigma: \Sigma estimate by PACE
 % along with 95% credible intervals
 % Btau, BT: basis function evaluations on tau and T.
 
%%    
addpath(genpath(cat(2, pwd, '/bspline')))  %include bspline package

t = sort(unique(cell2mat(T)));
% form the working grid
P = cellfun(@length, Y); % length of all observation grids
n = size(Y, 2); % # of signals
J = ones(p, 1);
D = abs(J * tau - tau' * J'); % distance matrix on the working grid tau
p_cgrid = length(t_cgrid);


%% MCMC sampling set up
%delta = 5
dm = delta + p - 1; %degrees of freedom for the inv-Wishart with delta
dmp = n + delta + p; %degrees of freedom for the conditional posterior inv-Wishart with n+1+delta
% c = 1; %Choose prior normal distribution for mean (arbitrarily)


%% Covariance estimate with PACE
if (isempty(Sigma_tau) || isempty(mu_tau_est))
    tau = t(floor(linspace(1, length(t), p))); % work grid
    addpath(genpath(cat(2, pwd, '/PACErelease2.11'))) %PACE package

    param_X = setOptions('regular', reg, 'FVE_threshold',1,'corrPlot',0,'rho',-1); 
    PCA_out= FPCA(Y, T, param_X);   %perform PCA on x, uncommon grid case

    Phihat=getVal(PCA_out,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(PCA_out,'lambda'); % The estimated eigenvalues.
    out1 = getVal(PCA_out, 'out1');
    Phi_tau = interp1(out1', Phihat, tau, 'spline');
    Sigma_tau = Phi_tau*diag(lamhat)*Phi_tau'; 

    mu_pca = getVal(PCA_out, 'mu');
    mu_tau_est = interp1(out1', mu_pca, tau, 'spline')';
else display('Run MCMC with given mean/covariance estimate on tau ... ')
end

Sigma_tau = topdm(Sigma_tau);

%% Solve for rho and nu
% by minimizing the MSE between A and Matern correlation matrix
if (mat) 
    
    A = COR(Sigma_tau);
    
    if(isempty(nu))       
        myobj = @(x) mean(mean((A - Matern(D, x(1), x(2), 1)).^2)); 
        [xx, fval] = fmincon(myobj, [1; 3], [], [], [], [], [1e-3; 2.5], [Inf; 100]);
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
  display(['minimum L2 distance=', num2str(fval), '; rho = ', num2str(rho), '; nu = ', num2str(nu)]);

else
    
    A = Sigma_tau;
    display('using empirical estimated covariance structure in the IW scale matrix');
    
end


%%
optknots = optknt(tau, 4);
Btau = bspline_basismatrix(4, optknots, tau);

BT_cgrid = bspline_basismatrix(4, optknots, t_cgrid);

BT = cell(1, n);
for i = 1: n
    BT{i} = bspline_basismatrix(4, optknots, T{i});
end 



%%
% Empirical estimates of noise variance and \sigma_s^2

func = @(x) sum(diff(x) .^ 2);
snhat2 =  sum(cellfun(func, Y)) / (2 * sum(P - 1)); 
rn = 1 / snhat2; % Sample estimated noise precision
rs = (trace(Sigma_tau) - p * snhat2 ) / (trace(A) / (delta - 2)); 

display(['Initial Estimates: ', ' noise variance = ', num2str(snhat2), ...
    '; sigma_s^2 = ', num2str(rs)])

%Choose inv-gamma distribution for noisy variance w = 1; ws = 5;
b = 1 / w; a = b * rn;
bs = 1 / ws; as = bs * rs;

%Choose gamma distribution for 1/noise-variance and signal-variance
display([' a = ',num2str(a),'; b = ',num2str(w),'; as = ',num2str(as),'; bs = ',num2str(ws)])


%%
% dimension p by 1
mu_tau_est = reshape(mu_tau_est, p, 1);
mu0 = mu_tau_est; 
mu_tau = mu_tau_est;

rn = 1 / snhat2; % Sample estimated noise precision
rs = (trace(Sigma_tau) - p * snhat2 ) / (trace(A) / (delta - 2)); 

Ztau = repmat(mu_tau, 1, n);
Zt = cell(1, n);

% Sigma_tau = Matern(D, 0.5, 3.5, sf^2);
% Sigma_tau = Phi_tau*diag(lamhat)*Phi_tau'; 
iK_tau = Sigma_tau;
K_tau = pinv(Sigma_tau);  % Covariance matrix Sigma on working grid

efunc = @(x, y) sum((x - y).^2);

% mean and covariance of coefficients zeta
mu_zeta = Btau \ mu_tau_est; % p X 1 dimension
iK_zeta = (Btau \ iK_tau) / Btau'; % p by p dimension
K_zeta = pinv(iK_zeta); 

Zeta = repmat(mu_zeta, 1, n);
mu0_zeta = Btau \ mu0;
A_zeta = (Btau\A) / (Btau'); 

%% assign memory ahead  to store MCMC samples
mu_cgrid_Out = NaN(p_cgrid, M); % Mean vector over t_cgrid
iK_cgrid_Out = NaN(p_cgrid, p_cgrid, M); % Covariance matrix over t_cgrid
Zt_cgrid_Out = NaN(p_cgrid, n, M); % Signal mcmc samples over t_cgrid

mu_zeta_Out = NaN(p, M); % Mean vector of zeta 
iK_zeta_Out = NaN(p, p, M); % Covariance matrix of zeta
ZetaOut = NaN(p, n, M); % Zeta mcmc samples
rnOut = NaN(1, M); % Precision of error, 1/noise-variance
rsOut = NaN(1, M); % signal variance, \sigma_s^2

%% Gibbs sampler
display('Start MCMC') 

for iter = 1 : (M + Burnin)    

    % update zeta/Z
        for i = 1:n
           % SB = iK_zeta * BT{i}';
           % BSB = BT{i} * SB;
           % SBI = SB /(BSB + (1/rn) * eye(P(i)));
           % mu_zeta_i = mu_zeta + SBI * (Y{i}' - BT{i} * mu_zeta);
           % iK_zeta_i = iK_zeta * Btau' - SBI * SB'; 
           
           iK_zeta_i = pinv(BT{i}' * BT{i} * rn + K_zeta);
           mu_zeta_i = iK_zeta_i * (BT{i}' * Y{i}' * rn + K_zeta * mu_zeta);
            
            Zeta(:, i) = mu_zeta_i + mychol(iK_zeta_i) * normrnd(0, 1, p, 1);
            Zt{i} = (BT{i} * Zeta(:, i))'; % approximate Z_i(ti)
        end
        
    %update 1/(noise variance), rn   
    rn = gamrnd(sum(P)/2 + a, 1/(b + sum(cellfun(efunc, Y, Zt)) / 2) );
    
    %update mu_zeta
    Mu_var = iK_zeta ./ (n + c);
    Mu_mean = (sum(Zeta, 2) + c .* mu0_zeta) ./ (n + c);    
    mu_zeta = Mu_mean + mychol(Mu_var) * normrnd(0, 1, p, 1);
    
    %update precision matrix K_zeta
    G = ((Zeta - repmat(mu_zeta, 1, n)) * (Zeta - repmat(mu_zeta, 1, n))' +...
        c .* (mu_zeta -  mu0_zeta) * (mu_zeta -  mu0_zeta)' + rs .* A_zeta);
    [iK_zeta, K_zeta] = myiwishrnd(dmp, G, p);
    
     %update signal variance rs
     iK_tau = (Btau) * iK_zeta * (Btau') ;
     K_tau = pinv(iK_tau) ;
     rs = gamrnd(p * dm/2 + as, 1/(bs + trace(A * K_tau) / 2));  

    % Save all MCMC samples   
    if iter > Burnin
       mu_zeta_Out(:, (iter - Burnin)) = mu_zeta;
       iK_zeta_Out(:, :, (iter - Burnin)) = iK_zeta; 
       ZetaOut(:, :, (iter - Burnin)) = Zeta;
       rnOut(iter-Burnin) = rn;
       rsOut(iter - Burnin) = rs;
       
       mu_cgrid_Out(:, (iter - Burnin)) = BT_cgrid * mu_zeta;
       iK_cgrid_Out(:, :, (iter - Burnin)) = BT_cgrid * iK_zeta * BT_cgrid';
       Zt_cgrid_Out(:, :, (iter - Burnin)) = (BT_cgrid * Zeta);
       
    end
    
end

display('End of MCMC ... ')


%% MCMC diagnosis
 addpath(genpath(cat(2, pwd, '/mcmcdiag')))
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
Zt_cgrid = mean(Zt_cgrid_Out(:, :, :), 3);

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

Zt_cgrid_sort = sort(Zt_cgrid_Out, 3); 
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

Zt_cgrid_CL = Zt_cgrid_sort(:, :, q1);
Zt_cgrid_UL = Zt_cgrid_sort(:, :, q2);
iK_cgrid_CL = iK_cgrid_sort(:, :, q1);
iK_cgrid_UL = iK_cgrid_sort(:, :, q2);
mu_cgrid_CI = [mu_cgrid_sort(:, q1), mu_cgrid_sort(:, q2)];

rs_CI = [rs_sort(q1), rs_sort(q2)];
rn_CI = [rn_sort(q1), rn_sort(q2)];


%%

if(mat)
    output = struct('tau', tau, 'Zt', {Zt}, 'iK_tau', iK_tau, 'iK_zeta_SE', iK_zeta_SE, ...
    'mu_tau', mu_tau, 'mu_zeta', mu_zeta, 'mu_tau_est', mu_tau_est, 'rn', rn, 'rs', rs,...
    'rho', rho, 'nu', nu, 'Sigma_tau', Sigma_tau, 'Zeta_CL', Zeta_CL, ...
    'Zeta_UL', Zeta_UL, 'iK_zeta_CL', iK_zeta_CL, ...
    'iK_zeta_UL', iK_zeta_UL, 'mu_zeta_CI', mu_zeta_CI, ...
    'Btau', Btau, 'BT', {BT}, 'Zeta', Zeta, 'iK_zeta', iK_zeta, ...
    'Zt_cgrid', Zt_cgrid, 'Zt_cgrid_CL', Zt_cgrid_CL, 'Zt_cgrid_UL', Zt_cgrid_UL, ...
    'iK_cgrid', iK_cgrid, 'iK_cgrid_CL', iK_cgrid_CL, 'iK_cgrid_UL', iK_cgrid_UL,...
    'mu_cgrid', mu_cgrid, 'mu_cgrid_CI', mu_cgrid_CI, 't_cgrid', t_cgrid, ...
    'rs_CI', rs_CI, 'rn_CI', rn_CI, 'mat', mat, 'optknots', optknots);
else
    output = struct('tau', tau, 'Zt', {Zt}, 'iK_tau', iK_tau, 'iK_zeta_SE', iK_zeta_SE, ...
    'mu_tau', mu_tau, 'mu_zeta', mu_zeta, 'mu_tau_est', mu_tau_est, 'rn', rn, 'rs', rs,...
    'Sigma_tau', Sigma_tau, 'Zeta_CL', Zeta_CL, ...
    'Zeta_UL', Zeta_UL, 'iK_zeta_CL', iK_zeta_CL, ...
    'iK_zeta_UL', iK_zeta_UL, 'mu_zeta_CI', mu_zeta_CI, ...
    'Btau', Btau, 'BT', {BT}, 'Zeta', Zeta, 'iK_zeta', iK_zeta, ...
    'Zt_cgrid', Zt_cgrid, 'Zt_cgrid_CL', Zt_cgrid_CL, 'Zt_cgrid_UL', Zt_cgrid_UL, ...
    'iK_cgrid', iK_cgrid, 'iK_cgrid_CL', iK_cgrid_CL, 'iK_cgrid_UL', iK_cgrid_UL,...
    'mu_cgrid', mu_cgrid, 'mu_cgrid_CI', mu_cgrid_CI, 't_cgrid', t_cgrid, ...
    'rs_CI', rs_CI, 'rn_CI', rn_CI, 'mat', mat, 'optknots', optknots);
end

display('Job completed.');




