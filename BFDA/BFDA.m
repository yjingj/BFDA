%% Author: Jingjing Yang (yjingj@gmail.com)

%% Main function to call BHM_mcmc, BABF_mcmc
% Bayesian GP regression will be added later...
%======
%Input: 
%======
%      Y:          1*n cell array, Y{i} is the vector of measurements for the ith subject,
%                  i=1,...,n.
%      T:          1*n cell array, T{i} is the vector of dependent variable (e.g. time points) 
%                    for the ith subject on which
%                  corresponding measurements Y{i} are taken, i=1,...,n.
%      param:      a struct obtained from setOptions.m sets the rest of arguments for BFDA.m
%                  e.g. 
%                  >> param = setOptions_bfda();
%=======
%Output:  
%=======  
%      out_smooth:         a data struct that contains all returned values
%                           of smoothing and mean-covariance estimation 
%   To see the names for out_smooth, type names(out_smooth)
%   To see an example, check with example_bfda.m
%
%%
function [ out_smooth, param ] = BFDA(Y, T, param)

  if nargin == 2
    param = setOptions_bfda(); % set default parameter options
  end
  
    pgrid = sort(unique(cell2mat(T))); %pooled grid
    p = length(pgrid);
    param.m = min(p, param.m);

    if strcmp(param.smethod, 'babf')
        if isempty(param.tau)
            param.tau = prctile(pgrid, 0:(100/(param.m-1)):100);
            % working grid
        end
    else
        param.tau = pgrid; % pooled grid
    end

%% Obtain pre-smooth estimator for mean and covariance on either working grid for babf or pooled grid for bhm 

if ~strcmp(param.smethod, 'css') && param.pace
    % By PACE 
    display('Run PACE to obtain pre-estimates for mean-covariance...')
    param_pace = setOptions('regular',param.cgrid,'FVE_threshold',0.9999,'corrPlot',0,'rho',-1);
    pace_out = FPCA(Y, T, param_pace); 
    
    Phihat=getVal(pace_out,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(pace_out,'lambda'); % The estimated eigenvalues.
    out1 = getVal(pace_out, 'out1');
    pace_mu = getVal(pace_out, 'mu');
    
    Phi_tau = interp1(out1', Phihat, param.tau, 'spline');
    param.Sigma_est = Phi_tau * diag(lamhat) * Phi_tau'; % Make sure positive definite 
    param.mu_est = interp1(out1', pace_mu, param.tau, 'spline');
    % Sigma_est = param.Sigma_est; mu_est = param.mu_est;
        
elseif ~strcmp(param.smethod, 'css') && (isempty(param.Sigma_est) || isempty(param.mu_est) )
    % Run bspline smoothing on raw data and then calculate the pooled
    % mean-covariance
    display('Run SplineSmoothing to obtain smoothed functional data for pre-estimates for mean-covariance...')
    n = length(Y);
    m = length(param.tau);
    Ysmooth = NaN(m, n);
    for i = 1:n
        Ysmooth(:, i) = csaps(T{i}, Y{i}, [], param.tau);
    end
    
    if isempty(param.Sigma_est)
        param.Sigma_est = cov(Ysmooth'); % Make sure positive definite
    end
    
    if isempty(param.mu_est)
        param.mu_est = mean(Ysmooth, 2);  
    end
    
elseif strcmp(param.smethod, 'css')
    display('CSS no need pre mean-covariance estimates.')
    
end

  %% smoothing and mean-covariance estimation
  
  if strcmp(param.smethod, 'babf')
      % with approximations by basis funcion
      display('Smoothing and estimation of functional data by BABF...')
      
      [out_smooth] = babf_mcmc(Y, T, param.delta, param.Burnin, param.M, param.mat, ...
          param.Sigma_est, param.mu_est, param.tau, ...
          param.w, param.ws, param.c, param.nu, param.eval_grid, param.resid_thin, param.tol);
      
  elseif strcmp(param.smethod, 'bhm')
      display('Smoothing and estimation of functional data by BHM...')
      [out_smooth] = bhm_mcmc(Y, T, param.delta, param.cgrid, param.Burnin, param.M, ...
                        param.mat, param.Sigma_est, param.mu_est, pgrid, param.nu, ...
                        param.c, param.w, param.ws, param.resid_thin, param.tol);
      
  elseif strcmp(param.smethod, 'bgp')
      display('Smoothing and estimation of functional data by BGP...')
    [out_smooth] = bgp(Y, T, param.Burnin, param.M, ...
                        param.mat, param.Sigma_est, param.mu_est,  ...
                        param.c, param.w, param.ws, param.nu, param.tol);
                    
  elseif strcmp(param.smethod, 'bfpca')
      display('Smoothing and estimation of functional data by BFPCA...')
    [out_smooth] = bfpca(Y, T, param.Burnin, param.M, ...
                        param.Sigma_est, param.a, ...
                        param.b, param.w);
                    
  elseif strcmp(param.smethod, 'css')
    % with BHM
    display('Smoothing and estimation of functional data by Cubic Smoothing Splines...')
    [out_smooth] = css_gcv(Y, T, param.lamb_min, param.lamb_max, param.lamb_step, param.eval_grid);
  else 
      display('param.smethod has to be one of: babf, bhm, bgp, bfpca, css. ')
    
  end
  

end

