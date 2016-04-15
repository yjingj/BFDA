%% Plot one sample functional curve
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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve'))

%% Plot another sample functional curve
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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve2'))

%% Plot CV errors
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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_fitted'))

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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_beta'))

%% plot predict responses with scalar responses

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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_pred'))

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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_fitted'))

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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_intercept'))

%%  plot the temperature coefficient function with functional responses
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell_fdy{2}, betastderrcell_fdy{2})
title('\beta(t) of Bayesian Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 2)
set(gca, 'fontsize', 14);
hold off

subplot(2, 1, 2)
plotbeta(betaestcell_raw_fdy{2}, betastderrcell_raw_fdy{2})
title('\beta(t) of Cubic Spline Smoothed Data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 2)
ylim([-1, 2.5])
set(gca, 'fontsize', 14);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_beta'))

%% plot predict responses with functional responses
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

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_pred'))

