%% Regression with smoothed functional data

% add fdaM path
addpath(genpath(cat(2, pwd, '/fdaM')))

% load Bayesian smoothing output
load('./Examples/Data/Simu_Data.mat')
load('./Examples/Data/Simu_Output.mat')

%%
n = 30; % Number of functional curves
p = 40; % Number of pooled grid points, or evaluated grid points
au = 0; bu = pi/2; % domain of t
pgrid = (au : (bu)/(p-1) : bu); % Pooled grid
trange = [au, bu];

sampind = sort(randsample(1:n,20,false)) ;
samptest = find(~ismember(1:n, sampind));
n_train = length(sampind); n_test = length(samptest);

cgrid = 0;

Xtrue = zeros(p, n);
Xraw = zeros(p, n);
Xsmooth = zeros(p, n);

if cgrid
% functional observations of the common-grid case 
    Xtrue = reshape(cell2mat(GausFD_cgrid.Xtrue_cell), p, n) ;
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

Xtrain = Xsmooth(:, sampind);
Xtest = Xsmooth(:, samptest);
Xraw_train = Xraw(:, sampind);
Xraw_test = Xraw(:, samptest);

rmse(Xtrue, Xsmooth)
rmse(Xtrue, Xraw)

%% Plot one sample curve
h = figure();
plot(GausFD_rgrid.Tcell{i}, GausFD_rgrid.Xraw_cell{i}, 'LineWidth', 3, 'Color', [0.75, 0.75, 0.75])
hold on
plot(pgrid, Xsmooth(:, i), 'r-',...
        pgrid, Xraw(:, i), 'b--', pgrid, Xtrue(:, i), 'c:', ...
        'LineWidth', 3, 'MarkerSize', 10)
xlabel('t'); ylabel('x(t)');
title('Example Functional Data')
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Raw Data', 'Bayesian Smoothed', 'Cubic Spline Smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve'))

%% Plot one sample curve
i=5;
h = figure();
plot(GausFD_rgrid.Tcell{i}, GausFD_rgrid.Xraw_cell{i}, 'LineWidth', 3, 'Color', [0.75, 0.75, 0.75])
hold on
plot(pgrid, Xsmooth(:, i), 'r-',...
        pgrid, Xraw(:, i), 'b--', pgrid, Xtrue(:, i), 'c:', ...
        'LineWidth', 3, 'MarkerSize', 10)
xlabel('t'); ylabel('x(t)');
title('Example Functional Data')
set(gca, 'fontsize', 14);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve2'))

%% response variables 
betamat = (pgrid') .^ 2 ;

% scalar respones
deltat  = pgrid(2)-pgrid(1);
Avec_true = deltat.*(Xtrue'*betamat    - ...
          0.5.*(Xtrue(1,    :)'.*betamat(1) + Xtrue(p,:)'.*betamat(p)) );
Avec = Avec_true + normrnd(0, 1, n, 1);

Avec_train = Avec(sampind);
Avec_test = Avec(samptest);
Avec_train_true = Avec_true(sampind);
Avec_test_true = Avec_true(samptest);

% functional responses
ymat_true = Xtrue .* repmat(betamat, 1, n) ;
ymat = ymat_true + normrnd(0, 1, p, n);

ymat_train = ymat(:, sampind);
ymat_test = ymat(:, samptest);
ymat_train_true = ymat_true(:, sampind);
ymat_test_true = ymat_true(:, samptest);

%%  set up xfd, yfd 
xnbasis = 20;
xbasis = create_bspline_basis(trange, xnbasis, 4);

xfd_true = smooth_basis(pgrid, Xtrue, xbasis);
xfd = smooth_basis(pgrid, Xtrain, xbasis);
xfd_raw = smooth_basis(pgrid, Xraw_train, xbasis);
[yfd_samp, df, gcv, beta, SSE, penmat, y2cMap, argvals, y] = smooth_basis(pgrid, ymat_train, xbasis);

%%  ----  set up the curvetrure penalty operator  -------
conbasis = create_constant_basis(trange); %  create a constant basis
wfd = fd([0, 1], conbasis);
wfdcell = fd2cell(wfd);
curvLfd = Lfd(2, wfdcell);

% set up xfdcell, betacell for scalar responses
xfdcell_vecy = cell(1, 2);
xfdcell_vecy{1} = fd(ones(1, n_train), conbasis);
xfdcell_vecy{2} = xfd;

xfd_raw_cell_vecy = cell(1, 2);
xfd_raw_cell_vecy{1} = fd(ones(1, n_train), conbasis);
xfd_raw_cell_vecy{2} = xfd_raw;

betafd0 = fd(0, conbasis);
bnbasis = 10;
betabasis = create_bspline_basis(trange, bnbasis, 4);
betafd1 = fd(zeros(bnbasis, 1), betabasis);

betacell_vecy = cell(1, 2);
betacell_vecy{1} = fdPar(betafd0);
betacell_vecy{2} = fdPar(betafd1, curvLfd, 0);

% set up xfdcell, betacell, yfd_par for functional responses
xfdcell_fdy = cell(1, 2);
xfdcell_fdy{1} = fd(ones(xnbasis, n_train), xbasis);
xfdcell_fdy{2} = xfd;

xfd_raw_cell_fdy = cell(1, 2);
xfd_raw_cell_fdy{1} = fd(ones(xnbasis, n_train), xbasis);
xfd_raw_cell_fdy{2} = xfd_raw;

betacell_fdy = cell(1, 2);
betacell_fdy{1} = fdPar(betafd1, curvLfd, 0);
betacell_fdy{2} = fdPar(betafd1, curvLfd, 0);

yfd_par = fdPar(yfd_samp, curvLfd, 0);

%%  compute cross-validated SSE's for a range of smoothing parameters
% %{ 
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
    SSE_CV_vecy(ilam) = fRegress_CV(Avec_train, xfdcell_vecy, betacelli_vecy, wt);
    fprintf('Scalar responses, lambda %6.2f: SSE = %10.4f\n', lam(ilam), SSE_CV_vecy(ilam));
    
    SSE_CV_raw_vecy(ilam) = fRegress_CV(Avec_train, xfd_raw_cell_vecy, betacelli_vecy, wt);
    fprintf('Scalar responses, lambda %6.2f: SSE = %10.4f\n', lam(ilam), SSE_CV_raw_vecy(ilam));
    
    betacelli_fdy    = betacell_fdy;
    betacelli_fdy{1} = putlambda(betacell_fdy{1}, lambda_vecy);
    betacelli_fdy{2} = putlambda(betacell_fdy{2}, lambda_vecy);
    yfd_par_i = putlambda(yfd_par, lambda_vecy);
    
    SSE_CV_fdy(ilam) = fRegress_CV(yfd_par_i, xfdcell_fdy, betacelli_fdy, wt);
    fprintf('Functional respones, lambda %6.2f: SSE = %10.4f\n', lam(ilam), SSE_CV_fdy(ilam));
    
    SSE_CV_raw_fdy(ilam) = fRegress_CV(yfd_par_i, xfd_raw_cell_fdy, betacelli_fdy, wt);
    fprintf('Functional respones, lambda %6.2f: SSE = %10.4f\n', lam(ilam), SSE_CV_raw_fdy(ilam));
end

figure()
subplot(2, 1, 1)
plot(lam, SSE_CV_vecy, 'ro-', lam, SSE_CV_raw_vecy, 'bo-')
xlabel('\fontsize{19} log_{10} smoothing parameter \lambda')
ylabel('\fontsize{19} Cross-validation score')
legend('Bayesian Smoothed', 'Cubic Spline Smoothed')
subplot(2, 1, 2)
plot(lam, SSE_CV_fdy, 'bo-', lam, SSE_CV_raw_fdy, 'ro-')
xlabel('\fontsize{19} log_{10} smoothing parameter \lambda')
ylabel('\fontsize{19} Cross-validation score')
% %}

%% Fit the linear model 
lambda = 0.1;
wt = ones(1, length(sampind));

% --------- scalar responses
betacell_vecy{2} = fdPar(betafd1, curvLfd, lambda);

fRegressStruct_vecy = fRegress(Avec_train, xfdcell_vecy, betacell_vecy, wt);
fRegressStruct_raw_vecy = fRegress(Avec_train, xfd_raw_cell_vecy, betacell_vecy, wt);

% get coefficients
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

% compute Rsquare
covmat = cov([Avec_train, Avec_hat]);
Rsqrd = covmat(1,2)^2/(covmat(1,1)*covmat(2,2));
disp(['R-squared = ',num2str(Rsqrd)])

covmat_raw = cov([Avec_train, Avec_hat_raw]);
Rsqrd_raw = covmat_raw(1,2)^2/(covmat_raw(1,1)*covmat_raw(2,2));
disp(['raw R-squared = ',num2str(Rsqrd_raw)])

%  compute sigma
resid_vecy = Avec_train - Avec_hat;
SigmaE_vecy = mean(resid_vecy.^2);
disp(['Scalar responses: SigmaE = ',num2str(SigmaE_vecy)])

resid_raw_vecy = Avec_train - Avec_hat_raw;
SigmaE_raw_vecy = mean(resid_raw_vecy.^2);
disp(['Scalar responses: Raw SigmaE = ',num2str(SigmaE_raw_vecy)])

% ---------- functional responses
betacell_fdy{1} = fdPar(betafd1, curvLfd, lambda);
betacell_fdy{2} = fdPar(betafd1, curvLfd, lambda);
yfd_par = fdPar(yfd_samp, curvLfd, lambda);

fRegressStruct_fdy = fRegress(yfd_par, xfdcell_fdy, betacell_fdy, wt, y2cMap);
fRegressStruct_raw_fdy = fRegress(yfd_par, xfd_raw_cell_fdy, betacell_fdy, wt, y2cMap);

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
 
%  compute squared residual correlation
resid_fdy = ymat_train_true - ymat_fitted;
SigmaE_fdy = cov(resid_fdy');

resid_raw_fdy = ymat_train_true - ymat_fitted_raw;
SigmaE_raw_fdy = cov(resid_raw_fdy');


%%  plot the fit with scalar repsonses

h= figure();
plot(Avec_train_true, Avec_hat, 'rx', ...
        Avec_train_true, Avec_hat_raw, 'bo', ...
        Avec_train_true, Avec_train_true, 'c--', 'MarkerSize', 10, 'Linewidth', 2)
ylabel(''); xlabel('True Area');
title('Fitted Response Values')
set(gca, 'fontsize', 16);
[hleg, hobj] = legend('Bayesian Smoothed', 'Cubic Spline Smoothed', 'Truth', 'Location', 'NorthWest');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 2); set(tobj, 'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_fitted'))

%% Plot fitted output with functional responses
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


%%  recompute the analysis to get confidence limits

% ------- scalar reponses
stderrStruct_vecy = fRegress_stderr(fRegressStruct_vecy, eye(n_train), SigmaE_vecy);
betastderrcell_vecy = stderrStruct_vecy.betastderr;

stderrStruct_raw_vecy = fRegress_stderr(fRegressStruct_raw_vecy, eye(n_train), SigmaE_raw_vecy);
betastderrcell_raw_vecy = stderrStruct_raw_vecy.betastderr;

%  constant  coefficient standard error:
intercept_vecy = getcoef(betastderrcell_vecy{1});
intercept_raw_vecy = getcoef(betastderrcell_raw_vecy{1});

% -------- functional responses
stderrStruct_fdy = fRegress_stderr(fRegressStruct_fdy, y2cMap, SigmaE_fdy);
% fixed a bug in fRegress_stderr.m at line 124: 
% bstderrfdj = data2fd(bstderrj, tfine, betabasisj); should be 
% bstderrfdj = data2fd(tfine, bstderrj, betabasisj);
betastderrcell_fdy = stderrStruct_fdy.betastderr;

stderrStruct_raw_fdy = fRegress_stderr(fRegressStruct_raw_fdy, y2cMap, SigmaE_raw_fdy);
betastderrcell_raw_fdy = stderrStruct_raw_fdy.betastderr;

%  coefficient standard error:
intercept_fdy = getcoef(betastderrcell_fdy{1});
intercept_raw_fdy = getcoef(betastderrcell_raw_fdy{1});


%%  plot the temperature coefficient function with scalar responses
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell_vecy{2}, betastderrcell_vecy{2})
ylim([-2.5, 3])
title('\beta(t) of Bayesian Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
set(gca, 'fontsize', 14);
hold off

subplot(2, 1, 2)
plotbeta(betaestcell_raw_vecy{2}, betastderrcell_raw_vecy{2})
ylim([-2.5, 3])
title('\beta(t) of Cubic Spline Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
set(gca, 'fontsize', 14);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_beta'))


%% plot intercept function with functional responses
h = figure();
plot(pgrid, intercept_fdy, 'r-', pgrid, intercept_raw_fdy, 'b--', 'LineWidth', 3);
xlabel('t');
ylabel('');
title('Intercept Function');
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Bayesian Smoothed', 'Cubic Spline Smoothed', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_intercept'))

%%  plot the temperature coefficient function with functional responses
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

Avec_pred = fRegressPred(xfdcell_test, betaestcell_vecy);

Avec_pred_raw = fRegressPred(xfd_raw_test_cell, betaestcell_raw_vecy);

display(['predict mse = ', num2str(mse(Avec_test_true, Avec_pred)), ...
    '; predict mse_raw = ',num2str(mse(Avec_test_true, Avec_pred_raw))])

%% plot predict responses

h=figure();
plot(Avec_test_true, Avec_pred, 'rx', ...
        Avec_test_true, Avec_pred_raw, 'bo', ...
        Avec_test_true, Avec_test_true, 'c--', 'MarkerSize', 10, 'Linewidth', 2)
ylabel(''); xlabel('Truth')
title('Predicted Response Values')
set(gca, 'fontsize', 16);
[hleg, hobj] = legend('Bayesian Smoothed', 'Cubic Spline Smoothed', 'Truth',...
                        'Location', 'NorthWest');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 2); set(tobj, 'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_pred'))

