% Regression with smoothed functional data by BFDA

%% Add fdaM path and load the functional data
addpath(genpath(cat(2, pwd, '/fdaM'))) % replace pwd by the directory of your MATLAB package
addpath(genpath(cat(2, pwd, '/BFDA_SourceCode')))

load('./Examples/Data/Simu_Data.mat')
load('./Examples/Data/Simu_Output.mat')

%% Set sample sizes, training data set, and test data set
n = 30; % Number of functional curves
p = 40; % Number of pooled grid points, or evaluated grid points
au = 0; bu = pi/2; % Domain of t
pgrid = (au : (bu)/(p-1) : bu); % Pooled grid
trange = [au, bu];

% sampind = sort(randsample(1:n,20,false)) ;
sampind = [1 2 3 4 5 6 7 9 11 12 13 15 17 19 22 23 24 25 28 30];
samptest = find(~ismember(1:n, sampind));
n_train = length(sampind); n_test = length(samptest);

cgrid = 0;

Xtrue = zeros(p, n);
Xraw = zeros(p, n);
Xsmooth = zeros(p, n);

if cgrid
% Functional observations with common grids 
    Xtrue = reshape(cell2mat(GausFD_cgrid.Xtrue_cell), p, n);
    Xsmooth = out_cgrid.Z(:, sampind);
    Xraw = reshape(cell2mat(GausFD_cgrid.Xraw_cell), p, n);
else
% Functional observations with random grids
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

Xtrain = Xsmooth(:, sampind);
Xtest = Xsmooth(:, samptest);
Xraw_train = Xraw(:, sampind);
Xraw_test = Xraw(:, samptest);

rmse(Xtrue, Xsmooth)
rmse(Xtrue, Xraw)

%% Generate response variables  
betamat = (pgrid') .^ 2 ;

% Scalar respones
deltat  = pgrid(2)-pgrid(1);
Avec_true = deltat.*(Xtrue'*betamat    - ...
          0.5.*(Xtrue(1,    :)'.*betamat(1) + Xtrue(p,:)'.*betamat(p)) );
Avec = Avec_true + normrnd(0, 1, n, 1);

Avec_train = Avec(sampind);
Avec_test = Avec(samptest);
Avec_train_true = Avec_true(sampind);
Avec_test_true = Avec_true(samptest);

% Functional responses
ymat_true = Xtrue .* repmat(betamat, 1, n) ;
ymat = ymat_true + normrnd(0, 1, p, n);

ymat_train = ymat(:, sampind);
ymat_test = ymat(:, samptest);
ymat_train_true = ymat_true(:, sampind);
ymat_test_true = ymat_true(:, samptest);

%%  Set up functional data structure xfd, yfd for fdaM
xnbasis = 20;
xbasis = create_bspline_basis(trange, xnbasis, 4);

xfd_true = smooth_basis(pgrid, Xtrue, xbasis);
xfd = smooth_basis(pgrid, Xtrain, xbasis);
xfd_raw = smooth_basis(pgrid, Xraw_train, xbasis);
[yfd_samp, df, gcv, beta, SSE, penmat, y2cMap, argvals, y] = ...
    smooth_basis(pgrid, ymat_train, xbasis);

%%  ----  Set up the curvatrure penalty operator  -------
conbasis = create_constant_basis(trange); %  create a constant basis
wfd = fd([0, 1], conbasis);
wfdcell = fd2cell(wfd);
curvLfd = Lfd(2, wfdcell);

% Set up xfdcell
xfdcell = cell(1, 2);
xfdcell{1} = fd(ones(1, n_train), conbasis);
xfdcell{2} = xfd;

xfd_raw_cell = cell(1, 2);
xfd_raw_cell{1} = fd(ones(1, n_train), conbasis);
xfd_raw_cell{2} = xfd_raw;

% Set up betacell for scalar responses
betafd0 = fd(0, conbasis);
bnbasis = 10;
betabasis = create_bspline_basis(trange, bnbasis, 4);
betafd1 = fd(zeros(bnbasis, 1), betabasis);

betacell_vecy = cell(1, 2);
betacell_vecy{1} = fdPar(betafd0);
betacell_vecy{2} = fdPar(betafd1, curvLfd, 0);

% Set up betacell, yfd_par for functional responses
betacell_fdy = cell(1, 2);
betacell_fdy{1} = fdPar(betafd1, curvLfd, 0);
betacell_fdy{2} = fdPar(betafd1, curvLfd, 0);

yfd_par = fdPar(yfd_samp, curvLfd, 0);

%%  Compute cross-validated SSE's for a range of smoothing parameters 
%{
wt = ones(1, length(sampind));
lam = (0:0.1:1);
nlam   = length(lam);

SSE_CV_vecy = zeros(nlam,1);
SSE_CV_raw_vecy = zeros(nlam, 1);

SSE_CV_fdy = zeros(nlam,1);
SSE_CV_raw_fdy = zeros(nlam, 1);

for ilam = 1:nlam;
   lambda_vecy       = lam(ilam);
    betacelli_vecy    = betacell_vecy;
    betacelli_vecy{2} = putlambda(betacell_vecy{2}, lambda_vecy);
    SSE_CV_vecy(ilam) = fRegress_CV(Avec_train, xfdcell, betacelli_vecy, wt);
    fprintf('Scalar responses, lambda %6.2f: SSE = %10.4f\n', ...
                lam(ilam), SSE_CV_vecy(ilam));
    
    SSE_CV_raw_vecy(ilam) = fRegress_CV(Avec_train, xfd_raw_cell, betacelli_vecy, wt);
    fprintf('Scalar responses, lambda %6.2f: SSE = %10.4f\n', ...
                lam(ilam), SSE_CV_raw_vecy(ilam));
    
    betacelli_fdy    = betacell_fdy;
    betacelli_fdy{1} = putlambda(betacell_fdy{1}, lambda_vecy);
    betacelli_fdy{2} = putlambda(betacell_fdy{2}, lambda_vecy);
    yfd_par_i = putlambda(yfd_par, lambda_vecy);
    
    SSE_CV_fdy(ilam) = fRegress_CV(yfd_par_i, xfdcell, betacelli_fdy, wt);
    fprintf('Functional respones, lambda %6.2f: SSE = %10.4f\n', ...
                lam(ilam), SSE_CV_fdy(ilam));
    
    SSE_CV_raw_fdy(ilam) = fRegress_CV(yfd_par_i, xfd_raw_cell, betacelli_fdy, wt);
    fprintf('Functional respones, lambda %6.2f: SSE = %10.4f\n', ...
                lam(ilam), SSE_CV_raw_fdy(ilam));
end
%}

%% Fit the linear model 
lambda = 0.1;
wt = ones(1, length(sampind));

% --------- Scalar responses
betacell_vecy{2} = fdPar(betafd1, curvLfd, lambda);

fRegressStruct_vecy = fRegress(Avec_train, xfdcell, betacell_vecy, wt);
fRegressStruct_raw_vecy = ...
    fRegress(Avec_train, xfd_raw_cell, betacell_vecy, wt);

% Get coefficients
betaestcell_vecy   = fRegressStruct_vecy.betahat; 
Avec_hat = fRegressStruct_vecy.yhat;
intercept_vecy = getcoef(getfd(betaestcell_vecy{1}));
disp(['Constant term = ',num2str(intercept_vecy)])

betaestcell_raw_vecy   = fRegressStruct_raw_vecy.betahat; 
Avec_hat_raw = fRegressStruct_raw_vecy.yhat;
intercept_raw_vecy = getcoef(getfd(betaestcell_raw_vecy{1}));
disp(['Constant term = ',num2str(intercept_raw_vecy)])

display(['Scalar reponses:', 'fitted mse = ', ...
    num2str(mse(Avec_train_true, Avec_hat)), ...
    '; fitted mse_raw = ',num2str(mse(Avec_train_true, Avec_hat_raw))])

% Compute Rsquare
covmat = cov([Avec_train, Avec_hat]);
Rsqrd = covmat(1,2)^2/(covmat(1,1)*covmat(2,2));
disp(['R-squared = ',num2str(Rsqrd)])

covmat_raw = cov([Avec_train, Avec_hat_raw]);
Rsqrd_raw = covmat_raw(1,2)^2/(covmat_raw(1,1)*covmat_raw(2,2));
disp(['raw R-squared = ',num2str(Rsqrd_raw)])

%  Compute sigma
resid_vecy = Avec_train - Avec_hat;
SigmaE_vecy = mean(resid_vecy.^2);
disp(['Scalar responses: SigmaE = ',num2str(SigmaE_vecy)])

resid_raw_vecy = Avec_train - Avec_hat_raw;
SigmaE_raw_vecy = mean(resid_raw_vecy.^2);
disp(['Scalar responses: Raw SigmaE = ',num2str(SigmaE_raw_vecy)])

% ---------- Functional responses
betacell_fdy{1} = fdPar(betafd1, curvLfd, lambda);
betacell_fdy{2} = fdPar(betafd1, curvLfd, lambda);
yfd_par = fdPar(yfd_samp, curvLfd, lambda);

fRegressStruct_fdy = fRegress(yfd_par, xfdcell, betacell_fdy, wt, y2cMap);
fRegressStruct_raw_fdy = ...
    fRegress(yfd_par, xfd_raw_cell, betacell_fdy, wt, y2cMap);

betaestcell_fdy   = fRegressStruct_fdy.betahat; 
yfd_hat = fRegressStruct_fdy.yhat;
intercept_fdy = eval_fd(pgrid, getfd(betaestcell_fdy{1}));

betaestcell_raw_fdy   = fRegressStruct_raw_fdy.betahat; 
yfd_hat_raw = fRegressStruct_raw_fdy.yhat;
intercept_raw_fdy = eval_fd(pgrid, getfd(betaestcell_raw_fdy{1}));

%  MSE of fitted responses 
ymat_fitted = eval_fd(pgrid, yfd_hat);
ymat_fitted_raw = eval_fd(pgrid, yfd_hat_raw);

display(['mse = ', num2str(mse(ymat_train_true, ymat_fitted)), ...
    '; mse_raw = ',num2str(mse(ymat_train_true, ymat_fitted_raw))])
 
%  Compute squared residual correlation
resid_fdy = ymat_train_true - ymat_fitted;
SigmaE_fdy = cov(resid_fdy');

resid_raw_fdy = ymat_train_true - ymat_fitted_raw;
SigmaE_raw_fdy = cov(resid_raw_fdy');

%%  recompute the analysis to get confidence limits
% ------- Scalar responses
stderrStruct_vecy = fRegress_stderr(fRegressStruct_vecy, eye(n_train), SigmaE_vecy);
betastderrcell_vecy = stderrStruct_vecy.betastderr;

stderrStruct_raw_vecy = ...
    fRegress_stderr(fRegressStruct_raw_vecy, eye(n_train), SigmaE_raw_vecy);
betastderrcell_raw_vecy = stderrStruct_raw_vecy.betastderr;

%  Constant  coefficient standard error:
intercept_std_vecy = getcoef(betastderrcell_vecy{1});
intercept_ste_raw_vecy = getcoef(betastderrcell_raw_vecy{1});

% -------- Functional responses
stderrStruct_fdy = fRegress_stderr(fRegressStruct_fdy, y2cMap, SigmaE_fdy);
% fixed a bug in fRegress_stderr.m at line 124: 
% bstderrfdj = data2fd(bstderrj, tfine, betabasisj); should be 
% bstderrfdj = data2fd(tfine, bstderrj, betabasisj);
betastderrcell_fdy = stderrStruct_fdy.betastderr;

stderrStruct_raw_fdy = ...
    fRegress_stderr(fRegressStruct_raw_fdy, y2cMap, SigmaE_raw_fdy);
betastderrcell_raw_fdy = stderrStruct_raw_fdy.betastderr;

%  Coefficient standard error:
intercept_std_fdy = eval_fd(pgrid, betastderrcell_fdy{1});
intercept_std_raw_fdy = eval_fd(pgrid, betastderrcell_raw_fdy{1});


%%  predict test data 
%  Set up xfd for test data
xfd_test = smooth_basis(pgrid, Xtest, xbasis);
xfd_raw_test = smooth_basis(pgrid, Xraw_test, xbasis);

% -------- Scalar responses
xfdcell_test = cell(1, 2);
xfdcell_test{1} = fd(ones(1, n_test), conbasis);
xfdcell_test{2} = xfd_test;

xfd_raw_test_cell = cell(1, 2);
xfd_raw_test_cell{1} = fd(ones(1, n_test), conbasis);
xfd_raw_test_cell{2} = xfd_raw_test;

Avec_pred = fRegressPred(xfdcell_test, betaestcell_vecy);
Avec_pred_raw = fRegressPred(xfd_raw_test_cell, betaestcell_raw_vecy);

display(['Scalar responses predict mse = ', ...
    num2str(mse(Avec_test_true, Avec_pred)), ...
    '; Scalar responses with raw data predict mse_raw = ',...
    num2str(mse(Avec_test_true, Avec_pred_raw))])

% -------- Functional responses
ymat_test_pred = ...
    eval_fd(pgrid, fRegressPred(xfdcell_test, betaestcell_fdy, xbasis));
ymat_test_pred_raw = ...
    eval_fd(pgrid, fRegressPred(xfd_raw_test_cell, betaestcell_raw_fdy, xbasis));

display(['Functional responses prediction mse = ', ...
    num2str(mse(ymat_test_true, ymat_test_pred)), ...
    '; Functional responses prediction with Raw data mse_raw = ',...
    num2str(mse(ymat_test_true, ymat_test_pred_raw))])

%% Use plot_script_fdregress.m to make plots
%close all;
%run('./Examples/plot_script_fdregress.m');
