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

if cgrid
% functional observations of the common-grid case 
    Xtrue = reshape(cell2mat(GausFD_cgrid.X), p, n) ;
    Xsmooth = out_cgrid.Z(:, sampind);
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

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve'))

%%  Call functional regression functions from fdaM library

% set up xfd 
xnbasis = 20;
xbasis = create_bspline_basis(trange, xnbasis, 4);

xfd_true = smooth_basis(pgrid, Xtrue, xbasis);
xfd = smooth_basis(pgrid, Xtrain, xbasis);
xfd_raw = smooth_basis(pgrid, Xraw_train, xbasis);

% response variable
nfine = p; tfine = pgrid';

deltat  = tfine(2)-tfine(1);
xmat    = eval_fd(tfine, xfd_true);
betamat = tfine .^ 2 ;
Avec_true = deltat.*(xmat'*betamat    - ...
          0.5.*(xmat(1,    :)'.*betamat(1) + xmat(p,:)'.*betamat(p)) );
Avec = Avec_true + normrnd(0, 1, n, 1);

Avec_train = Avec(sampind);
Avec_test = Avec(samptest);

Avec_train_true = Avec_true(sampind);
Avec_test_true = Avec_true(samptest);


%%  ----  set up the harmonic acceleration, curvetrure penalty operator  -------
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
xfdcell{1} = fd(ones(1, n_train), conbasis);
xfdcell{2} = xfd;

xfd_raw_cell = cell(1, 2);
xfd_raw_cell{1} = fd(ones(1, n_train), conbasis);
xfd_raw_cell{2} = xfd_raw;

% set up betacell
betafd0 = fd(0, conbasis);

bnbasis = 10;
% betabasis = create_fourier_basis(trange, bnbasis);
betabasis = create_bspline_basis(trange, bnbasis, 4);
betafd1 = fd(zeros(bnbasis, 1), betabasis);

betacell = cell(1, 2);
betacell{1} = fdPar(betafd0);
% betacell{2} = fdPar(betafd1, harmaccelLfd, 0);
betacell{2} = fdPar(betafd1, curvLfd, 0);

%%  compute cross-validated SSE's for a range of smoothing parameters
%{ 
wt = ones(1, length(sampind));
lam = (0:0.1:1);
nlam   = length(lam);
SSE_CV = zeros(nlam,1);
SSE_CV_raw = zeros(nlam, 1);
for ilam = 1:nlam;
   lambda       = lam(ilam);
    betacelli    = betacell;
    betacelli{2} = putlambda(betacell{2}, lambda);
    SSE_CV(ilam) = fRegress_CV(Avec_train, xfdcell, betacelli, wt);
    SSE_CV_raw(ilam) = fRegress_CV(Avec_train, xfd_raw_cell, betacelli, wt);
    fprintf('%3.f %6.2f %10.4f\n', ilam, lam(ilam), SSE_CV(ilam));
    fprintf('%3.f %6.2f %10.4f\n', ilam, lam(ilam), SSE_CV_raw(ilam));
end

plot(lam, SSE_CV, 'ro-', lam, SSE_CV_raw, 'bo-')
xlabel('\fontsize{19} log_{10} smoothing parameter \lambda')
ylabel('\fontsize{19} Cross-validation score')
legend('Smoothed data', 'Raw data')
%}

%% Fit the linear model 

lambda = 0.1;
betacell{2} = fdPar(betafd1, curvLfd, lambda);

% call regression function 
wt = ones(1, length(sampind));
fRegressStruct = fRegress(Avec_train, xfdcell, betacell, wt);
fRegressStruct_raw = fRegress(Avec_train, xfd_raw_cell, betacell, wt);

betaestcell   = fRegressStruct.betahat; 
Avec_hat = fRegressStruct.yhat;
intercept = getcoef(getfd(betaestcell{1}));
disp(['Constant term = ',num2str(intercept)])

betaestcell_raw   = fRegressStruct_raw.betahat; 
Avec_hat_raw = fRegressStruct_raw.yhat;
intercept_raw = getcoef(getfd(betaestcell_raw{1}));
disp(['Constant term = ',num2str(intercept_raw)])

display(['fitted mse = ', num2str(mse(Avec_train_true, Avec_hat)), ...
    '; fitted mse_raw = ',num2str(mse(Avec_train_true, Avec_hat_raw))])

%%  plot the fit

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


%%  compute squared multiple correlation

covmat = cov([Avec_train, Avec_hat]);
Rsqrd = covmat(1,2)^2/(covmat(1,1)*covmat(2,2));
disp(['R-squared = ',num2str(Rsqrd)])

covmat_raw = cov([Avec_train, Avec_hat_raw]);
Rsqrd_raw = covmat_raw(1,2)^2/(covmat_raw(1,1)*covmat_raw(2,2));
disp(['raw R-squared = ',num2str(Rsqrd_raw)])

%  compute sigma
resid = Avec_train - Avec_hat;
SigmaE = mean(resid.^2);
disp(['SigmaE = ',num2str(SigmaE)])

resid_raw = Avec_train - Avec_hat_raw;
SigmaE_raw = mean(resid_raw.^2);
disp(['raw SigmaE = ',num2str(SigmaE_raw)])

%%  recompute the analysis to get confidence limits

%  get confidence intervals

stderrStruct = fRegress_stderr(fRegressStruct, eye(n_train), SigmaE);

betastderrcell = stderrStruct.betastderr;

stderrStruct_raw = fRegress_stderr(fRegressStruct_raw, eye(n_train), SigmaE_raw);

betastderrcell_raw = stderrStruct_raw.betastderr;

%  constant  coefficient standard error:

intercept = getcoef(betastderrcell{1});
intercept_raw = getcoef(betastderrcell_raw{1});

%%  plot the temperature coefficient function
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell{2}, betastderrcell{2})
ylim([-2.5, 3])
title('\beta(t) of Bayesian Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 2)
set(gca, 'fontsize', 16);
hold off

subplot(2, 1, 2)
plotbeta(betaestcell_raw{2}, betastderrcell_raw{2})
ylim([-2.5, 3])
title('\beta(t) of Cubic Spline Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 2)
set(gca, 'fontsize', 16);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_beta'))

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

Avec_pred = fRegressPred(xfdcell_test, betaestcell);

Avec_pred_raw = fRegressPred(xfd_raw_test_cell, betaestcell_raw);

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

