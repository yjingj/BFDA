% Example plot script for the BFDA simulation and smoothing

%% Load BFDA smoothed data
% close all

load('./Examples/Data/Simu_Data.mat')
load('./Examples/Data/Simu_Output.mat')

%%
[xgrid,ygrid] = meshgrid(eval_grid); % 3D plots of the covariance matrices

clims = [0, 1];
xlims = [0, pi/2];

%% set up data for stationary-cgrid-ucgrid plot

Xraw_cgrid = reshape(cell2mat(GausFD_cgrid.Xraw_cell), length(GausFD_cgrid.Xraw_cell{1}), n);
Xtrue_cgrid = reshape(cell2mat(GausFD_cgrid.Xtrue_cell), length(GausFD_cgrid.Xtrue_cell{1}), n);
T_cgrid = reshape(cell2mat(GausFD_cgrid.Tcell), length(GausFD_cgrid.Tcell{1}), n);

Xfull_cgrid = Xraw_cgrid;
mu_sample_cgrid = nanmean(Xfull_cgrid, 2);


Xraw_ucgrid = reshape(cell2mat(GausFD_ucgrid.Xraw_cell), length(GausFD_ucgrid.Xraw_cell{1}), n);
Xtrue_ucgrid = reshape(cell2mat(GausFD_ucgrid.Xtrue_cell), length(GausFD_ucgrid.Xtrue_cell{1}), n);
T_ucgrid = reshape(cell2mat(GausFD_ucgrid.Tcell), length(GausFD_ucgrid.Tcell{1}), n);

Xfull_ucgrid = NaN(p, n); % n by p data matrix with nan's for unobserved data
    for i = 1:n
       Idx = find(ismember(pgrid, GausFD_ucgrid.Tcell{i}));
       Xfull_ucgrid(Idx, i) = GausFD_ucgrid.Xraw_cell{i};
    end
mu_sample_ucgrid = nanmean(Xfull_ucgrid, 2);

%% plot stationary functional data with common and uncommon grids

leglocation='northeast';
ylims = [-10, 8]; 

h = figure('DefaultLegendFontSize',18,'DefaultLegendFontSizeMode','manual');

subplot_tight(2, 2, 1,  [.1, .05])

idx = [1];
plot(T_cgrid(:, idx), Xraw_cgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_cgrid(:, idx), Xtrue_cgrid(:, idx), 'r-.', eval_grid, out_cgrid.Z(:, idx), 'b-', eval_grid, out_cgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

[hleg1, hobj1] = legend('Sample', 'Truth', 'BHM', 'BHM 95% CI', 'Location', leglocation);

idx = [6];
plot(T_cgrid(:, idx), Xraw_cgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_cgrid(:, idx), Xtrue_cgrid(:, idx), 'r-.', eval_grid, out_cgrid.Z(:, idx), 'b-', eval_grid, out_cgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(a)'); %Functional data (common grids).
set(gca,'fontsize', 20)
%htext = findobj(get(hleg1, 'children'), 'type', 'text');
%set(htext,'FontSize', 38)
%legend boxoff


subplot_tight(2, 2, 2,  [.1, .05])

idx = [1];
plot(T_ucgrid(:, idx), Xraw_ucgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_ucgrid(:, idx), Xtrue_ucgrid(:, idx), 'r-.', eval_grid, out_ucgrid.Z(:, idx), 'b-', eval_grid, out_ucgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

idx = [6];
plot(T_ucgrid(:, idx), Xraw_ucgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_ucgrid(:, idx), Xtrue_ucgrid(:, idx), 'r-.', eval_grid, out_ucgrid.Z(:, idx), 'b-', eval_grid, out_ucgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(b)');
%title('(b) Functional data (uncommon grids).');
set(gca,'fontsize', 20)

subplot_tight(2, 2, 3,  [.1, .05])
plot(eval_grid, mu_sample_cgrid, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_cgrid.Mean_true, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_cgrid.mu,  'b-', eval_grid, out_cgrid.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(c)');
%title('(c) Means (common grids).');
xlim(xlims)
ylim([-5, 5])
set(gca,'fontsize', 20)

subplot_tight(2, 2, 4,  [.1, .05])
plot(eval_grid, mu_sample_ucgrid, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_ucgrid.Mean_true, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_ucgrid.mu,  'b-', eval_grid, out_ucgrid.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(d)');
%title('(d) Means (uncommon grids).');
xlim(xlims)
ylim([-5, 5])
set(gca,'fontsize', 20)


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/stat_sig'))

%% plot covariance 

zlims = [-1, 8];

h = figure();

subplot_tight(2, 2, 1,  [.1,.05])
mesh(xgrid, ygrid, out_cgrid.Sigma)
zlim(zlims)
title('(a)'); % BHM covariance estimate(common grids).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 2,  [.1,.05])
mesh(xgrid, ygrid, out_ucgrid.Sigma)
zlim(zlims)
title('(b)'); %BHM covariance estimate(uncommon grids).
set(gca,'fontsize', 20)


subplot_tight(2, 2, 3,  [.1,.05])
mesh(xgrid, ygrid, cov(Xfull_cgrid'))
zlim(zlims)
title('(c)'); %Sample covariance estimate(common grids).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 4,  [.1,.05])
mesh(xgrid, ygrid,  GausFD_ucgrid.Cov_true)
zlim(zlims)
title('(d)'); %True covariance(common/uncommon grids).
set(gca,'fontsize', 20)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/stat_cov'))


%% set up data for Nonstationary-cgrid-ucgrid plot

Xraw_cgrid_ns = reshape(cell2mat(GausFD_cgrid_ns.Xraw_cell), length(GausFD_cgrid_ns.Xraw_cell{1}), n);
Xtrue_cgrid_ns = reshape(cell2mat(GausFD_cgrid_ns.Xtrue_cell), length(GausFD_cgrid_ns.Xtrue_cell{1}), n);
T_cgrid_ns = reshape(cell2mat(GausFD_cgrid_ns.Tcell), length(GausFD_cgrid_ns.Tcell{1}), n);

Xfull_cgrid_ns = Xraw_cgrid_ns;
mu_sample_cgrid_ns = nanmean(Xfull_cgrid_ns, 2);


Xraw_ucgrid_ns = reshape(cell2mat(GausFD_ucgrid_ns.Xraw_cell), length(GausFD_ucgrid_ns.Xraw_cell{1}), n);
Xtrue_ucgrid_ns = reshape(cell2mat(GausFD_ucgrid_ns.Xtrue_cell), length(GausFD_ucgrid_ns.Xtrue_cell{1}), n);
T_ucgrid_ns = reshape(cell2mat(GausFD_ucgrid_ns.Tcell), length(GausFD_ucgrid_ns.Tcell{1}), n);

Xfull_ucgrid_ns = NaN(p, n); % n by p data matrix with nan's for unobserved data
    for i = 1:n
       Idx = find(ismember(pgrid, GausFD_ucgrid_ns.Tcell{i}));
       Xfull_ucgrid_ns(Idx, i) = GausFD_ucgrid_ns.Xraw_cell{i};
    end
mu_sample_ucgrid_ns = nanmean(Xfull_ucgrid_ns, 2);

%% plot stationary functional data with common and uncommon grids

leglocation='southwest';
ylims = [-18, 18]; 

h = figure('DefaultLegendFontSize',18,'DefaultLegendFontSizeMode','manual');

subplot_tight(2, 2, 1,  [.1, .05])

idx = [1];
plot(T_cgrid_ns(:, idx), Xraw_cgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_cgrid_ns(:, idx), Xtrue_cgrid_ns(:, idx), 'r-.', eval_grid, out_cgrid_ns.Z(:, idx), 'b-', ...
    eval_grid, out_cgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

[hleg1, hobj1] = legend('Sample', 'Truth', 'BHM', 'BHM 95% CI', 'Location', leglocation);

idx = [6];
plot(T_cgrid_ns(:, idx), Xraw_cgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_cgrid_ns(:, idx), Xtrue_cgrid_ns(:, idx), 'r-.', eval_grid, out_cgrid_ns.Z(:, idx), 'b-', ...
    eval_grid, out_cgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(a)'); % Functional data (common grids).
set(gca,'fontsize', 20)
%set(hleg1,'FontSize',20, 'color','none')
%legend boxoff


subplot_tight(2, 2, 2,  [.1, .05])

idx = [1];
plot(T_ucgrid_ns(:, idx), Xraw_ucgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_ucgrid_ns(:, idx), Xtrue_ucgrid_ns(:, idx), 'r-.', eval_grid, out_ucgrid_ns.Z(:, idx), 'b-',...
    eval_grid, out_ucgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

idx = [10];
plot(T_ucgrid_ns(:, idx), Xraw_ucgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_ucgrid_ns(:, idx), Xtrue_ucgrid_ns(:, idx), 'r-.', eval_grid, out_ucgrid_ns.Z(:, idx), 'b-', ...
    eval_grid, out_ucgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(b)'); % Functional data (uncommon grids).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 3,  [.1, .05])
plot(eval_grid, mu_sample_cgrid_ns, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_cgrid_ns.Mean_true, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_cgrid_ns.mu,  'b-', eval_grid, out_cgrid_ns.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(c)'); % Means (common grids).
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 20)

subplot_tight(2, 2, 4,  [.1, .05])
plot(eval_grid, mu_sample_ucgrid_ns, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_ucgrid_ns.Mean_true, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_ucgrid_ns.mu,  'b-', eval_grid, out_ucgrid_ns.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(d)'); % Means (uncommon grids).
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 20)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/non_stat_sig'))

%%
zlims = [-1, 35];

h = figure();

subplot_tight(2, 2, 1,  [.1,.05])
mesh(xgrid, ygrid, out_cgrid_ns.Sigma)
zlim(zlims)
title('(a)'); % BHM covariance estimate (common grids).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 2,  [.1,.05])
mesh(xgrid, ygrid, out_ucgrid_ns.Sigma)
zlim(zlims)
title('(b)');
set(gca,'fontsize', 20) % BHM covariance estimate (uncommon grids).


subplot_tight(2, 2, 3,  [.1,.05])
mesh(xgrid, ygrid, cov(Xfull_cgrid_ns'))
zlim(zlims)
title('(c)');
set(gca,'fontsize', 20) % Sample covariance estimate (common grids).

subplot_tight(2, 2, 4,  [.1,.05])
mesh(xgrid, ygrid,  GausFD_ucgrid_ns.Cov_true)
zlim(zlims)
title('(d)');
set(gca,'fontsize', 20) % True covariance (common/uncommon grids).

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/non_stat_cov'))


%% Setup data for random grids case

Xraw_rgrid = reshape(cell2mat(GausFD_rgrid.Xraw_cell), length(GausFD_rgrid.Xraw_cell{1}), n);
Xtrue_rgrid = reshape(cell2mat(GausFD_rgrid.Xtrue_cell), length(GausFD_rgrid.Xtrue_cell{1}), n);
T_rgrid = reshape(cell2mat(GausFD_rgrid.Tcell), length(GausFD_rgrid.Tcell{1}), n);

Xraw_rgrid_ns = reshape(cell2mat(GausFD_rgrid_ns.Xraw_cell), length(GausFD_rgrid_ns.Xraw_cell{1}), n);
Xtrue_rgrid_ns = reshape(cell2mat(GausFD_rgrid_ns.Xtrue_cell), length(GausFD_rgrid_ns.Xtrue_cell{1}), n);
T_rgrid_ns = reshape(cell2mat(GausFD_rgrid_ns.Tcell), length(GausFD_rgrid_ns.Tcell{1}), n);


%% plot stationary/nonstationary functional data with random grids

leglocation='southwest';
ylims = [-11, 11]; 

h = figure('DefaultLegendFontSize',18,'DefaultLegendFontSizeMode','manual');

subplot_tight(2, 2, 1,  [.1,.05])

idx = [1];
plot(T_rgrid(:, idx), Xraw_rgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_rgrid(:, idx), Xtrue_rgrid(:, idx), 'r-.', eval_grid, out_rgrid.Z_cgrid(:, idx), 'b-', ...
    eval_grid, out_rgrid.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

[hleg1, hobj1] = legend('Sample', 'Truth', 'BABF', 'BABF 95% CI', 'Location', leglocation);

idx = [6];
plot(T_rgrid(:, idx), Xraw_rgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_rgrid(:, idx), Xtrue_rgrid(:, idx), 'r-.', eval_grid, out_rgrid.Z_cgrid(:, idx), 'b-', ...
    eval_grid, out_rgrid.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(a)'); % Functional data (stationary).
set(gca,'fontsize', 20)
%set(hleg1,'FontSize',12, 'color','none')
%legend boxoff


subplot_tight(2, 2, 2,  [.1,.05])

idx = [1];
plot(T_rgrid_ns(:, idx), Xraw_rgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_rgrid_ns(:, idx), Xtrue_rgrid_ns(:, idx), 'r-.', eval_grid, out_rgrid_ns.Z_cgrid(:, idx), 'b-',...
    eval_grid, out_rgrid_ns.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid_ns.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

idx = [10];
plot(T_rgrid_ns(:, idx), Xraw_rgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_rgrid_ns(:, idx), Xtrue_rgrid_ns(:, idx), 'r-.', eval_grid, out_rgrid_ns.Z_cgrid(:, idx), 'b-', ...
    eval_grid, out_rgrid_ns.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid_ns.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim([-10, 13])
title('(b)'); % Functional data (nonstationary).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 3,  [.1,.05])
plot(GausFD_rgrid.Tcell{1}, GausFD_rgrid.Mean_true_cell{1}, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_rgrid.mu_cgrid,  'b-', eval_grid, out_rgrid.mu_cgrid_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(c)'); % Mean estimate (stationary).
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 20)

subplot_tight(2, 2, 4,  [.1,.05])
plot(GausFD_rgrid_ns.Tcell{1}, GausFD_rgrid_ns.Mean_true_cell{1}, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_rgrid_ns.mu_cgrid,  'b-', eval_grid, out_rgrid_ns.mu_cgrid_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(d)'); % Mean estimate (nonstationary).
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 20)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/rand_sig'))

%%
zlims = [-1, 25];

h = figure();

subplot_tight(2, 2, 1,  [.1,.05])
mesh(xgrid, ygrid, out_rgrid.Sigma_cgrid)
zlim([-1, 8])
title('(a)'); % BHM covariance estimate (stationary).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 2,  [.1,.05])
mesh(xgrid, ygrid, out_rgrid_ns.Sigma_cgrid)
zlim(zlims)
title('(b)'); % BHM covariance estimate (nonstationary).
set(gca,'fontsize', 20)


subplot_tight(2, 2, 3,  [.1,.05])
mesh(xgrid, ygrid, GausFD_cgrid.Cov_true)
zlim([-1, 8])
title('(c)'); % True covariance (stationary).
set(gca,'fontsize', 20)

subplot_tight(2, 2, 4,  [.1,.05])
mesh(xgrid, ygrid,  GausFD_cgrid_ns.Cov_true)
zlim(zlims)
title('(d)'); % True covariance (nonstationary).
set(gca,'fontsize', 20)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

%print(h, '-dpdf', cat(2, getenv('HOME'), '/Figures/rand_cov'))





