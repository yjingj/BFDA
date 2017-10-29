% Example plot script for the functional regression results

%%
xlims = [0, pi/2];

%% Plot one sample functional curve
i=30;
h = figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
subplot_tight(1, 1, 1,  [.1, .1])

plot(GausFD_rgrid.Tcell{i}, GausFD_rgrid.Xraw_cell{i}, 'LineWidth', 3, 'Color', [0.75, 0.75, 0.75])
hold on
plot(pgrid, Xsmooth(:, i), 'r-',...
        pgrid, Xraw(:, i), 'b--', pgrid, Xtrue(:, i), 'c:', ...
        'LineWidth', 3, 'MarkerSize', 10)
xlabel('t'); ylabel('x(t)');
title('(a)') %Example functional data
xlim(xlims)
[hleg, hobj] = legend('Raw data', 'BABF', 'CSS', 'Truth', 'Location', 'Best');
set(gca, 'fontsize', 26);

hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_sample_curve'))


%%  plot the fit with scalar repsonses

h= figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
subplot_tight(1, 1, 1,  [.1, .1])

plot(Avec_train_true, Avec_hat, 'rx', ...
        Avec_train_true, Avec_hat_raw, 'bo', ...
        Avec_train_true, Avec_train_true, 'c-', 'MarkerSize', 10, 'Linewidth', 3)
ylabel(''); xlabel('Truth'); ylabel('Fitted');
title('(a)') % Fitted response values
xlim([-5, 2]);
[hleg, hobj] = legend('BABF', 'CSS', 'Location', 'NorthWest');

set(gca, 'fontsize', 26);

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_vecy_fitted'))

%%  plot the temperature coefficient function with scalar responses
h = figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
subplot_tight(2, 1, 1,  [.1, .07])
plotbeta(betaestcell_vecy{2}, betastderrcell_vecy{2})
ylim([-2.5, 3])
title('(b)') %\beta(t) of Bayesian smoothed data
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
set(gca, 'fontsize', 26);
hold off

subplot_tight(2, 1, 2,  [.1, .07])
plotbeta(betaestcell_raw_vecy{2}, betastderrcell_raw_vecy{2})
ylim([-2.5, 3])
title('(c)') % \beta(t) of cubic spline smoothed data
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
set(gca, 'fontsize', 26);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_vecy_beta'))

%% plot predict responses with scalar responses

h=figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
%h=figure();
subplot_tight(1, 1, 1,  [.1, .1])

plot(Avec_test_true, Avec_pred, 'rx', ...
        Avec_test_true, Avec_pred_raw, 'bo', ...
        Avec_test_true, Avec_test_true, 'c-', 'MarkerSize', 10, 'Linewidth', 3)
ylabel(''); xlabel('Truth'); ylabel('Predicted');
title('(b)') %Predicted response values
xlim([-5, 3]);
[hleg, hobj] = legend('BABF', 'CSS', 'Location', 'NorthWest');

set(gca, 'fontsize', 26);

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_vecy_pred'))

%% Plot fitted output with functional responses
h=figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
i=10;

subplot_tight(1, 1, 1,  [.1, .1])
plot(pgrid, ymat_fitted(:, i), 'r-', ...
        pgrid, ymat_fitted_raw(:, i), 'b--', ...
        pgrid, ymat_train_true(:, i), 'c:', 'LineWidth', 3)
    
xlabel('t'); ylabel(''); xlim(xlims);
title('(a)') % fitted functional responses
[hleg, hobj] = legend('BABF', 'CSS', 'Truth', 'Location', 'SouthWest');
set(gca, 'fontsize', 26);
% set(hleg, 'fontsize', 20);
% lobj = findobj(hobj, 'type', 'line');
% tobj = findobj(hobj, 'type', 'text');
% set(lobj, 'linewidth', 3); set(tobj, 'fontsize', 14)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_fdy_fitted'))

%% plot intercept function with functional responses
h = figure('DefaultLegendFontSize',18,'DefaultLegendFontSizeMode','manual');
subplot_tight(2, 1, 1,  [.1, .07])
plotbeta(betaestcell_fdy{1}, betastderrcell_fdy{1})
title('(a)') % Intercept \beta_0(t) of Bayesian smoothed data
ylim([-1, 1]);
set(gca, 'fontsize', 26);

subplot_tight(2, 1, 2,  [.1, .07])
plotbeta(betaestcell_raw_fdy{1}, betastderrcell_raw_fdy{1})
title('(c)') % Intercept \beta_0(t) of cubic spline smoothed data
ylim([-1, 1]);
set(gca, 'fontsize', 26);

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_fdy_intercept'))

%%  plot the temperature coefficient function with functional responses
h = figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
subplot_tight(2, 1, 1,  [.1, .07])
plotbeta(betaestcell_fdy{2}, betastderrcell_fdy{2})
title('(b)') %\beta(t) of Bayesian smoothed data
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth',3)
set(gca, 'fontsize', 26);
ylim([-0.2, 3]);
hold off

subplot_tight(2, 1, 2,  [.1, .07])
plotbeta(betaestcell_raw_fdy{2}, betastderrcell_raw_fdy{2})
title('(d)') % \beta(t) of cubic spline smoothed data
hold on 
plot(pgrid, pgrid.^2, 'c:', 'Linewidth', 3)
ylim([-0.2, 3]);
set(gca, 'fontsize', 26);
hold off

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_fdy_beta'))

%% plot predict responses with functional responses
h=figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
i=5;
subplot_tight(1, 1, 1,  [.1, .1])
plot(pgrid, ymat_test_pred(:, i), 'r-', ...
        pgrid, ymat_test_pred_raw(:, i), 'b--', ...
        pgrid, ymat_test_true(:, i), 'c:', 'LineWidth', 3)
    
xlabel('t'); ylabel(''); xlim(xlims);
title('(b)') % Predicted functional response
set(gca, 'fontsize', 26);
[hleg, hobj] = legend('BABF', 'CSS', 'Truth', 'Location', 'SouthWest');

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/reg_fdy_pred'))

