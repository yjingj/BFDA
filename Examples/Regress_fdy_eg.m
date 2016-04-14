%% Regression with smoothed functional data by BFDA

% add fdaM path
addpath(genpath(cat(2, pwd, '/fdaM')))

% load Bayesian smoothing output
load('./Examples/Data/Simu_Data.mat')
load('./Examples/Data/Simu_Output.mat')

%% Get functionl covariate

n = 30; % Number of functional curves
p = 40; % Number of pooled grid points, or evaluated grid points
au = 0; bu = pi/2; % domain of t
pgrid = (au : (bu)/(p-1) : bu); % Pooled grid, i.e., evalutation grid
trange = [au, bu];

cgrid = 0;
if cgrid
% functional observations of the common-grid case 
    Xtrue = reshape(cell2mat(GausFD_cgrid.X), p, n) ;
    Xsmooth = out_cgrid.Z;
    Xraw = reshape(cell2mat(GausFD_cgrid.Y), p, n);
else
% functional observations of the random-grid case
    Xtrue = zeros(p, n);
    Xraw = zeros(p, n);
    for i = 1:n
        xi = GausFD_rgrid.X{i};
        xrawi = GausFD_rgrid.Y{i};
        ti = GausFD_rgrid.T{i};
        %interpolate by cubic smoothing spline
        h = mean(diff(ti));
        Xtrue(:, i) = csaps(ti, xi, 1/(1 + h^3/6), pgrid); 
        Xraw(:, i) = csaps(ti, xrawi, 1/(1 + h^3/6), pgrid);
        zi = out_rgrid.Zt{i};
        Xsmooth(:, i) = csaps(ti, zi, 1/(1 + h^3/6), pgrid);
    end
end

sampind = sort(randsample(1:n,20,false)) ;
samptest = find(~ismember(1:n, sampind));
n_train = length(sampind); n_test = length(samptest);

Xtrain = Xsmooth(:, sampind);
Xtest = Xsmooth(:, samptest);
Xraw_train = Xraw(:, sampind);
Xraw_test = Xraw(:, samptest);

rmse(Xtrue, Xsmooth)
rmse(Xtrue, Xraw)

%% Plot one sample curve
i=3;
h = figure();
plot(GausFD_rgrid.T{i}, GausFD_rgrid.Y{i}, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75])
hold on
plot(pgrid, Xsmooth(:, i), 'r-',...
        pgrid, Xraw(:, i), 'b--', pgrid, Xtrue(:, i), 'c:', ...
        'LineWidth', 2, 'MarkerSize', 10)
xlabel('t'); ylabel('x(t)');
title('One Sample Curve')
set(gca, 'fontsize', 16);
[hleg, hobj] = legend('Raw Data', 'Bayesian Smoothed', 'Cubic Spline Smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 2); set(tobj, 'fontsize', 16)
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve2'))


%%  Generate functional responses
betamat = (pgrid') .^ 2 ;
ymat_true = Xtrue .* repmat(betamat, 1, n) ;
ymat = ymat_true + normrnd(0, 1, p, n);

ymat_train = ymat(:, sampind);
ymat_test = ymat(:, samptest);
ymat_train_true = ymat_true(:, sampind);
ymat_test_true = ymat_true(:, samptest);


%%  Call functional regression functions from fdaM library

% set up xfd, yfd

xnbasis = 20;
xbasis = create_bspline_basis(trange, xnbasis, 4); % same basis for x(t) and response(t)

xfd_true = smooth_basis(pgrid, Xtrue, xbasis);

xfd = smooth_basis(pgrid, Xtrain, xbasis);

xfd_raw = smooth_basis(pgrid, Xraw_train, xbasis);

[yfd_samp, df, gcv, beta, SSE, penmat, y2cMap, argvals, y] = smooth_basis(pgrid, ymat_train, xbasis);

%%  ----  set up the harmonic acceleration operator  -------
% Lcoef   = [0,(2*pi/period)^2,0]; 
conbasis = create_constant_basis(trange); %  create a constant basis
Lcoef   = [0,16,0];    %  set up three coefficients
wfd     = fd(Lcoef,conbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

wfd = fd([0, 1], conbasis);
wfdcell = fd2cell(wfd);
curvLfd = Lfd(2, wfdcell);

% set up xfdcell
xfdcell = cell(1, 2);
xfdcell{1} = fd(ones(xnbasis, n_train), xbasis);
xfdcell{2} = xfd;

xfd_raw_cell = cell(1, 2);
xfd_raw_cell{1} = fd(ones(xnbasis, n_train), xbasis);
xfd_raw_cell{2} = xfd_raw;

% set up betacell

bnbasis = 10;
% betabasis = create_fourier_basis(trange, bnbasis);
betabasis = create_bspline_basis(trange, bnbasis, 4);
betafd0 = fd(zeros(bnbasis, 1), betabasis);
betafd1 = fd(zeros(bnbasis, 1), betabasis);

betacell = cell(1, 2);
betacell{1} = fdPar(betafd0, curvLfd, 0);
% betacell{2} = fdPar(betafd1, harmaccelLfd, 0);
betacell{2} = fdPar(betafd1, curvLfd, 0);

yfd_par = fdPar(yfd_samp, curvLfd, 0);

%%  compute cross-validated SSE's for a range of smoothing parameters
wt = ones(1, length(sampind));
lam = (0:0.1:1);
nlam   = length(lam);
SSE_CV = zeros(nlam,1);
SSE_CV_raw = zeros(nlam, 1);
for ilam = 1:nlam;
    lambda       = lam(ilam);
    betacelli    = betacell;
    betacelli{1} = putlambda(betacell{1}, lambda);
    betacelli{2} = putlambda(betacell{2}, lambda);
    yfd_par_i = putlambda(yfd_par, lambda);
    
    SSE_CV(ilam) = fRegress_CV(yfd_par_i, xfdcell, betacelli, wt);
    fprintf('%3.f %6.2f %10.4f\n', ilam, lam(ilam), SSE_CV(ilam));
    
    SSE_CV_raw(ilam) = fRegress_CV(yfd_par_i, xfd_raw_cell, betacelli, wt);
    fprintf('%3.f %6.2f %10.4f\n', ilam, lam(ilam), SSE_CV_raw(ilam));
end

plot(lam, SSE_CV, 'bo-', lam, SSE_CV_raw, 'ro-')
xlabel('\fontsize{19} log_{10} smoothing parameter \lambda')
ylabel('\fontsize{19} Cross-validation score')

%% Fit the linear model 

lambda = 0.2;
betacell{1} = fdPar(betafd0, curvLfd, lambda);
betacell{2} = fdPar(betafd1, curvLfd, lambda);
yfd_par = fdPar(yfd_samp, curvLfd, lambda);

% call regression function 
wt = ones(1, n_train);
fRegressStruct = fRegress(yfd_par, xfdcell, betacell, wt, y2cMap);
fRegressStruct_raw = fRegress(yfd_par, xfd_raw_cell, betacell, wt, y2cMap);

%% fit output
betaestcell   = fRegressStruct.betahat; 
yfd_hat = fRegressStruct.yhat;
intercept = eval_fd(pgrid, getfd(betaestcell{1}));

betaestcell_raw   = fRegressStruct_raw.betahat; 
yfd_hat_raw = fRegressStruct_raw.yhat;
intercept_raw = eval_fd(pgrid, getfd(betaestcell_raw{1}));

%  MSE of fitted responses 
ymat_fitted = eval_fd(pgrid, yfd_hat);
ymat_fitted_raw = eval_fd(pgrid, yfd_hat_raw);

display(['mse = ', num2str(mse(ymat_train_true, ymat_fitted)), ...
    '; mse_raw = ',num2str(mse(ymat_train_true, ymat_fitted_raw))])
 
%  compute squared residual correlation
resid = ymat_train_true - ymat_fitted;
SigmaE = cov(resid');

resid_raw = ymat_train_true - ymat_fitted_raw;
SigmaE_raw = cov(resid_raw');

%% Plot fitted output
h=figure();
i=5;
plot(pgrid, ymat_fitted(:, i), 'r-', ...
        pgrid, ymat_fitted_raw(:, i), 'b--', ...
        pgrid, ymat_train_true(:, i), 'c:', 'MarkerSize', 16, 'LineWidth', 2)
    
xlabel('t'); ylabel('');
title('Fitted Functional Response')
set(gca, 'fontsize', 16);
[hleg, hobj] = legend('Bayesian Smoothed', 'Cubic Spline Smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 2); set(tobj, 'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_fitted'))


%% plot intercept function 
h = figure();
plot(pgrid, intercept, 'r-', pgrid, intercept_raw, 'b--', 'LineWidth', 2);
xlabel('t');
ylabel('');
title('Intercept Function');
set(gca, 'fontsize', 16);
[hleg, hobj] = legend('Bayesian Smoothed', 'Cubic Spline Smoothed', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 2); set(tobj, 'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_intercept'))


%%  recompute the analysis to get confidence limits

stderrStruct = fRegress_stderr(fRegressStruct, y2cMap, SigmaE);
% fixed a bug in fRegress_stderr.m at line 124: 
% bstderrfdj = data2fd(bstderrj, tfine, betabasisj); should be 
% bstderrfdj = data2fd(tfine, bstderrj, betabasisj);

betastderrcell = stderrStruct.betastderr;

stderrStruct_raw = fRegress_stderr(fRegressStruct_raw, y2cMap_raw, SigmaE_raw);

betastderrcell_raw = stderrStruct_raw.betastderr;

%  constant  coefficient standard error:

intercept = getcoef(betastderrcell{1});
intercept_raw = getcoef(betastderrcell_raw{1});

%%  plot the temperature coefficient function
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell{2}, betastderrcell{2})
title('\beta(t) of Bayesian Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 2)
set(gca, 'fontsize', 16);
hold off

subplot(2, 1, 2)
plotbeta(betaestcell_raw{2}, betastderrcell_raw{2})
title('\beta(t) of Cubic Spline Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 2)
ylim([-1, 2.5])
set(gca, 'fontsize', 16);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_beta'))


%%  predict test data 

%  set up functional data object for test data
xfd_test = smooth_basis(pgrid, Xtest, xbasis);

xfd_raw_test = smooth_basis(pgrid, Xraw_test, xbasis);

% set up xfdcell
xfdcell_test = cell(1, 2);
xfdcell_test{1} = fd(ones(1, n_test), conbasis);
xfdcell_test{2} = xfd_test;

xfd_raw_test_cell = cell(1, 2);
xfd_raw_test_cell{1} = fd(ones(1, n_test), conbasis);
xfd_raw_test_cell{2} = xfd_raw_test;

ymat_test_pred = eval_fd(pgrid, fRegressPred(xfdcell_test, betaestcell, xbasis));

ymat_test_pred_raw = eval_fd(pgrid, fRegressPred(xfd_raw_test_cell, betaestcell_raw, xbasis));

display(['mse = ', num2str(mse(ymat_test_true, ymat_test_pred)), ...
    '; mse_raw = ',num2str(mse(ymat_test_true, ymat_test_pred_raw))])

%% plot predict responses

h=figure();
i=2;
plot(pgrid, ymat_test_pred(:, i), 'r-', ...
        pgrid, ymat_test_pred_raw(:, i), 'b--', ...
        pgrid, ymat_test_true(:, i), 'c:', 'MarkerSize', 16, 'LineWidth', 2)
    
xlabel('t'); ylabel('');
title('Predicted Functional Response')
set(gca, 'fontsize', 16);
[hleg, hobj] = legend('Bayesian Smoothed', 'Cubic Spline Smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 2); set(tobj, 'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_pred'))



