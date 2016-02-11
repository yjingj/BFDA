% Main function to call BHM_mcmc, BABF_mcmc, Regress
%======
%Input: 
%======
%      Y:          1*n cell array, Y{i} is the vector of measurements for the ith subject,
%                  i=1,...,n.
%      T:          1*n cell array, T{i} is the vector of dependent variable (e.g. time points) 
%                    for the ith subject on which
%                  corresponding measurements Y{i} are taken, i=1,...,n.
%      param:          a struct obtained from setOptions.m sets the rest of arguments for BFDA.m
%                  e.g.
%                  >> param = setOptions_bfda();
%=======
%Output:  
%=======  
%      out_str:         a data struct that contains all returned values 
% 
%   To see the names for out_str, type names(out_str)
%   To see an example, check with example_bfda.m
%

function [ out_str ] = BFDA(Y, T, param)

  if nargin == 2
    param = setOptions_bfda(); % set default parameter options
  end
  
    pgrid = sort(unique(cell2mat(T))); %pooled grid

    if param.babf 
        tau = param.tau; % working grid
    else
        tau = pgrid; % pooled grid
    end
%% Obtain pre-smooth estimator for mean and covariance on either working grid for babf or pooled grid for bhm 

if param.pace
    % By PACE, add PACE path
    display('Run PACE to obtain pre-estimates for mean-covariance...')
    addpath(genpath(cat(2, pwd, '/PACErelease2.11')))
    param_pace = setOptions('regular',param.cgrid,'FVE_threshold',0.99,'corrPlot',0,'rho',-1);
    pace_out = FPCA(Y, T, param_pace); 
    
    Phihat=getVal(pace_out,'phi');  % The estimated Eigenfunctions.
    lamhat=getVal(pace_out,'lambda'); % The estimated eigenvalues.
    out1 = getVal(pace_out, 'out1');
    pace_mu = getVal(pace_out, 'mu');
    
    Phi_tau = interp1(out1', Phihat, tau, 'spline');
    param.Sigma_est = topdm(Phi_tau * diag(lamhat) * Phi_tau'); % Make sure positive definite 
    param.mu_est = interp1(out1', pace_mu, tau, 'spline');
else
    % Run bspline smoothing on raw data and then calculate the pooled
    % mean-covariance
    display('Run SplineSmoothing to obtain smoothed functional data for pre-estimates for mean-covariance...')
    n = length(Y);
    m = length(tau);
    Ysmooth = NaN(m, n);
    for i = 1:n
        Ysmooth(:, i) = csaps(T{i}, Y{i}, [], tau);
    end
    param.Sigma_est = cov(Ysmooth'); % Make sure positive definite 
    param.mu_est = mean(Ysmooth, 2);   
end

  %% smoothing and mean-covariance estimation
  if param.babf
      % with approximations by basis funcion
      display('Smoothing and estimation of functional data by BABF...')
      [out_str] = babf_mcmc(Y, T, param.delta, param.Burnin, param.M, param.mat, ...
          param.Sigma_est, param.mu_est, param.tau, param.p, ...
          param.w, param.ws, param.c, param.nu, param.reg, pgrid);
      
  else
    % with BHM
    display('Smoothing and estimation of functional data by BHM...')
    [out_str] = bhm_mcmc(Y, T, param.delta, param.cgrid, param.Burnin, param.M, ...
                        param.mat, param.Sigma_est, param.mu_est, pgrid, param.nu, ...
                        param.c, param.w, param.ws);
  end
  
  if param.Regress
      % regression 
      display('Regression with smoothed functional data ... ')
      
  end


end

