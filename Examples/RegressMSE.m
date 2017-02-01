% Replicate the regression example as in FdRegress_eg.m for 
% 100 times to calculate the mean square error of predictions.

%% set a random seed
stream = RandStream('twister','Seed', 1234);
reset(stream);  % set up a seed for simulation

%% add fdaM path and load Bayesian smoothing output
addpath(genpath(cat(2, pwd, '/fdaM'))) % replace pwd by the directory of your MATLAB package

load('./Examples/Data/Simu_Data.mat')
load('./Examples/Data/Simu_Output.mat')

%% setup data 
n = 30; % Number of functional curves
p = 40; % Number of pooled grid points, or evaluated grid points
au = 0; bu = pi/2; % domain of t
pgrid = (au : (bu)/(p-1) : bu); % Pooled grid
trange = [au, bu];

cgrid = 0;

Xtrue = zeros(p, n);
Xraw = zeros(p, n);
Xsmooth = zeros(p, n);

if cgrid
% functional observations of the common-grid case 
    Xtrue = reshape(cell2mat(GausFD_cgrid.Xtrue_cell), p, n);
    Xsmooth = out_cgrid.Z(:, sampind);
    Xraw = reshape(cell2mat(GausFD_cgrid.Xraw_cell), p, n);
else
% functional observations of the random-grid case
    for i = 1:n
        xi = GausFD_rgrid.Xtrue_cell{i};
        xrawi = GausFD_rgrid.Xraw_cell{i};
        ti = GausFD_rgrid.Tcell{i};
        %interpolate by cubic smoothing spline
        h = mean(diff(ti));
        Xtrue(:, i) = csaps(ti, xi, 1/(1 + h^3/6), pgrid); 
        Xraw(:, i) = csaps(ti, xrawi, 1/(1 + h^3/6), pgrid);
        zi = out_rgrid.Zt{i};
        Xsmooth(:, i) = csaps(ti, zi, 1/(1 + h^3/6), pgrid);
    end
    % Xsmooth = out_rgrid.Z_cgrid;
    % Xsmooth(1, :) = Xsmooth(2, :) * 0.9; 
end
rmse(Xtrue, Xsmooth)
rmse(Xtrue, Xraw)

betamat = (pgrid') .^ 2 ;

%% response variables 
% scalar respones
deltat  = pgrid(2)-pgrid(1);
Avec_true = deltat.*(Xtrue'*betamat    - ...
          0.5.*(Xtrue(1,    :)'.*betamat(1) + Xtrue(p,:)'.*betamat(p)) );
Avec = Avec_true + normrnd(0, 1, n, 1);

ymat_true = Xtrue .* repmat(betamat, 1, n) ;
ymat = ymat_true + normrnd(0, 1, p, n);

xnbasis = 20;
xbasis = create_bspline_basis(trange, xnbasis, 4);
xfd_true = smooth_basis(pgrid, Xtrue, xbasis);

%%  ----  set up the curvetrure penalty operator  -------
conbasis = create_constant_basis(trange); %  create a constant basis
wfd = fd([0, 1], conbasis);
wfdcell = fd2cell(wfd);
curvLfd = Lfd(2, wfdcell);

betafd0 = fd(0, conbasis);
bnbasis = 10;
betabasis = create_bspline_basis(trange, bnbasis, 4);
betafd1 = fd(zeros(bnbasis, 1), betabasis);

lambda = 0.1;
wt = ones(1, 20);

% set up betacell for scalar responses
betacell_vecy = cell(1, 2);
betacell_vecy{1} = fdPar(betafd0);
betacell_vecy{2} = fdPar(betafd1, curvLfd, lambda);

% set up betacell, yfd_par for functional responses
betacell_fdy = cell(1, 2);
betacell_fdy{1} = fdPar(betafd1, curvLfd, lambda);
betacell_fdy{2} = fdPar(betafd1, curvLfd, lambda);

%% Repeat the following study for 100 times
M = 100;

MSE_fit = zeros(4, M);
MSE_test = zeros(4, M);
SAMP_IND = zeros(M, 20);

%% 
for j = 1: M
    
sampind = sort(randsample(1:n,20,false)) ;
SAMP_IND(j, :) = sampind;
samptest = find(~ismember(1:n, sampind));
n_train = length(sampind); n_test = length(samptest);

Xtrain = Xsmooth(:, sampind);
Xtest = Xsmooth(:, samptest);
Xraw_train = Xraw(:, sampind);
Xraw_test = Xraw(:, samptest);

% response variables 
Avec_train = Avec(sampind);
Avec_test = Avec(samptest);
Avec_train_true = Avec_true(sampind);
Avec_test_true = Avec_true(samptest);

% functional responses
ymat_train = ymat(:, sampind);
ymat_test = ymat(:, samptest);
ymat_train_true = ymat_true(:, sampind);
ymat_test_true = ymat_true(:, samptest);

%  set up xfd, yfd 
xfd = smooth_basis(pgrid, Xtrain, xbasis);
xfd_raw = smooth_basis(pgrid, Xraw_train, xbasis);
[yfd_samp, df, gcv, beta, SSE, penmat, y2cMap, argvals, y] = smooth_basis(pgrid, ymat_train, xbasis);

% set up xfdcell
xfdcell = cell(1, 2);
xfdcell{1} = fd(ones(1, n_train), conbasis);
xfdcell{2} = xfd;

xfd_raw_cell = cell(1, 2);
xfd_raw_cell{1} = fd(ones(1, n_train), conbasis);
xfd_raw_cell{2} = xfd_raw;

yfd_par = fdPar(yfd_samp, curvLfd, 0);

% Fit the linear model 

% --------- scalar responses
fRegressStruct_vecy = fRegress(Avec_train, xfdcell, betacell_vecy, wt);
fRegressStruct_raw_vecy = fRegress(Avec_train, xfd_raw_cell, betacell_vecy, wt);

% get coefficients
betaestcell_vecy   = fRegressStruct_vecy.betahat; 
Avec_hat = fRegressStruct_vecy.yhat;

betaestcell_raw_vecy   = fRegressStruct_raw_vecy.betahat; 
Avec_hat_raw = fRegressStruct_raw_vecy.yhat;

MSE_fit(1, j) = mse(Avec_train_true, Avec_hat);
MSE_fit(2, j) = mse(Avec_train_true, Avec_hat_raw);

% ---------- functional responses
yfd_par = fdPar(yfd_samp, curvLfd, lambda);
fRegressStruct_fdy = fRegress(yfd_par, xfdcell, betacell_fdy, wt, y2cMap);
fRegressStruct_raw_fdy = fRegress(yfd_par, xfd_raw_cell, betacell_fdy, wt, y2cMap);

betaestcell_fdy   = fRegressStruct_fdy.betahat; 
yfd_hat = fRegressStruct_fdy.yhat;

betaestcell_raw_fdy   = fRegressStruct_raw_fdy.betahat; 
yfd_hat_raw = fRegressStruct_raw_fdy.yhat;

%  MSE of fitted responses 
ymat_fitted = eval_fd(pgrid, yfd_hat);
ymat_fitted_raw = eval_fd(pgrid, yfd_hat_raw);

MSE_fit(3, j) = mse(ymat_train_true, ymat_fitted);
MSE_fit(4, j) = mse(ymat_train_true, ymat_fitted_raw);
 
%  predict test data 
%  set up xfd for test data
xfd_test = smooth_basis(pgrid, Xtest, xbasis);
xfd_raw_test = smooth_basis(pgrid, Xraw_test, xbasis);

xfdcell_test = cell(1, 2);
xfdcell_test{1} = fd(ones(1, n_test), conbasis);
xfdcell_test{2} = xfd_test;

xfd_raw_test_cell = cell(1, 2);
xfd_raw_test_cell{1} = fd(ones(1, n_test), conbasis);
xfd_raw_test_cell{2} = xfd_raw_test;

% -------- scalar responses
Avec_pred = fRegressPred(xfdcell_test, betaestcell_vecy);
Avec_pred_raw = fRegressPred(xfd_raw_test_cell, betaestcell_raw_vecy);

MSE_test(1, j) = mse(Avec_test_true, Avec_pred);
MSE_test(2, j) = mse(Avec_test_true, Avec_pred_raw);

% -------- functional responses
ymat_test_pred = eval_fd(pgrid, fRegressPred(xfdcell_test, betaestcell_fdy, xbasis));
ymat_test_pred_raw = eval_fd(pgrid, fRegressPred(xfd_raw_test_cell, betaestcell_raw_fdy, xbasis));

MSE_test(3, j) = mse(ymat_test_true, ymat_test_pred);
MSE_test(4, j) = mse(ymat_test_true, ymat_test_pred_raw);

end

%% average MSE (mean square error) and std (standard deviation)
display('MSE of the fitted responses: Bayesian smoothed (scalar response), cubic spline smoothed (scalar response)')
display('1) Bayesian smoothed (scalar response)')
display('2) cubic spline smoothed (scalar response)')
display('3) Bayesian smoothed (functional response)')
display('4) cubic spline smoothed (functional response)')
mean(MSE_fit, 2)

display('std of the fitted responses')
std(MSE_fit, 0, 2)

display('MSE of the prediction responses')
mean(MSE_test, 2)
display('std of the prediction responses')
std(MSE_test, 0, 2)
