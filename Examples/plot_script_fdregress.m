%% Plot one sample functional curve
h = figure();
plot(GausFD_rgrid.Tcell{i}, GausFD_rgrid.Xraw_cell{i}, 'LineWidth', 3, 'Color', [0.75, 0.75, 0.75])
hold on
plot(pgrid, Xsmooth(:, i), 'r-',...
        pgrid, Xraw(:, i), 'b--', pgrid, Xtrue(:, i), 'c:', ...
        'LineWidth', 3, 'MarkerSize', 10)
xlabel('t'); ylabel('x(t)');
title('Example functional data')
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Raw data', 'Bayesian smoothed', 'Cubic spline smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

% print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve'))

%% Plot another sample functional curve
i=5;
h = figure();
plot(GausFD_rgrid.Tcell{i}, GausFD_rgrid.Xraw_cell{i}, 'LineWidth', 3, 'Color', [0.75, 0.75, 0.75])
hold on
plot(pgrid, Xsmooth(:, i), 'r-',...
        pgrid, Xraw(:, i), 'b--', pgrid, Xtrue(:, i), 'c:', ...
        'LineWidth', 3, 'MarkerSize', 10)
xlabel('t'); ylabel('x(t)');
title('Example functional data')
set(gca, 'fontsize', 14);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_sample_curve2'))


%%  plot the fit with scalar repsonses

h= figure();
plot(Avec_train_true, Avec_hat, 'rx', ...
        Avec_train_true, Avec_hat_raw, 'bo', ...
        Avec_train_true, Avec_train_true, 'c--', 'MarkerSize', 10, 'Linewidth', 3)
ylabel(''); xlabel('True Area');
title('Fitted response values')
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Bayesian smoothed', 'Cubic spline smoothed', 'Truth', 'Location', 'NorthWest');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_fitted'))

%%  plot the temperature coefficient function with scalar responses
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell_vecy{2}, betastderrcell_vecy{2})
ylim([-2.5, 3])
title('\beta(t) of Bayesian smoothed data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
set(gca, 'fontsize', 14);
hold off

subplot(2, 1, 2)
plotbeta(betaestcell_raw_vecy{2}, betastderrcell_raw_vecy{2})
ylim([-2.5, 3])
title('\beta(t) of cubic spline smoothed data')
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
        Avec_test_true, Avec_test_true, 'c--', 'MarkerSize', 10, 'Linewidth', 3)
ylabel(''); xlabel('Truth')
title('Predicted response values')
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Bayesian smoothed', 'Cubic spline smoothed', 'Truth',...
                        'Location', 'NorthWest');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_vecy_pred'))

%% Plot fitted output with functional responses
h=figure();
i=5;
plot(pgrid, ymat_fitted(:, i), 'r-', ...
        pgrid, ymat_fitted_raw(:, i), 'b--', ...
        pgrid, ymat_train_true(:, i), 'c:', 'LineWidth', 3)
    
xlabel('t'); ylabel('');
title('Fitted functional response')
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Bayesian smoothed', 'Cubic spline smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_fitted'))

%% plot intercept function with functional responses
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell_fdy{1}, betastderrcell_fdy{1})
title('Intercept \beta_0(t) of Bayesian smoothed data')
set(gca, 'fontsize', 14);
subplot(2, 1, 2)
plotbeta(betaestcell_raw_fdy{1}, betastderrcell_raw_fdy{1})
title('Intercept \beta_0(t) of cubic spline smoothed data')
set(gca, 'fontsize', 14);

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_intercept'))

%%  plot the temperature coefficient function with functional responses
h = figure();
subplot(2, 1, 1)
plotbeta(betaestcell_fdy{2}, betastderrcell_fdy{2})
title('\beta(t) of Bayesian smoothed data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth',3)
set(gca, 'fontsize', 14);
hold off

subplot(2, 1, 2)
plotbeta(betaestcell_raw_fdy{2}, betastderrcell_raw_fdy{2})
title('\beta(t) of cubic spline smoothed data')
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
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
        pgrid, ymat_test_true(:, i), 'c:', 'LineWidth', 3)
    
xlabel('t'); ylabel('');
title('Predicted functional response')
set(gca, 'fontsize', 14);
[hleg, hobj] = legend('Bayesian smoothed', 'Cubic spline smoothed', 'Truth', 'Location', 'Best');
set(hleg, 'fontsize', 20);
lobj = findobj(hobj, 'type', 'line');
tobj = findobj(hobj, 'type', 'text');
set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/reg_fdy_pred'))

