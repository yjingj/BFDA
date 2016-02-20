% close all
load('./Examples/Data/Simu_Data.mat')
load('./Examples/Data/Simu_Output.mat')

%%
n=30; p = 40; 
pgrid = (0 : (pi/2)/(p-1) : (pi/2));
eval_grid = (0 : (pi/2)/(p-1) : (pi/2));
[xgrid,ygrid] = meshgrid(eval_grid); % 3D plots of the covariance matrices

clims = [0, 1];
xlims = [0, pi/2];


%% set up data for stationary-cgrid-ucgrid plot

Y_cgrid = reshape(cell2mat(GausFD_cgrid.Y), length(GausFD_cgrid.Y{1}), n);
X_cgrid = reshape(cell2mat(GausFD_cgrid.X), length(GausFD_cgrid.X{1}), n);
T_cgrid = reshape(cell2mat(GausFD_cgrid.T), length(GausFD_cgrid.T{1}), n);

Yfull_cgrid = Y_cgrid;
mu_sample_cgrid = nanmean(Yfull_cgrid, 2);


Y_ucgrid = reshape(cell2mat(GausFD_ucgrid.Y), length(GausFD_ucgrid.Y{1}), n);
X_ucgrid = reshape(cell2mat(GausFD_ucgrid.X), length(GausFD_ucgrid.X{1}), n);
T_ucgrid = reshape(cell2mat(GausFD_ucgrid.T), length(GausFD_ucgrid.T{1}), n);

Yfull_ucgrid = NaN(p, n); % n by p data matrix with nan's for unobserved data
    for i = 1:n
       Idx = find(ismember(pgrid, GausFD_ucgrid.T{i}));
       Yfull_ucgrid(Idx, i) = GausFD_ucgrid.Y{i};
    end
mu_sample_ucgrid = nanmean(Yfull_ucgrid, 2);

%% plot stationary functional data with common and uncommon grids

leglocation='northeast';
ylims = [-8, 8]; 

h = figure();

subplot_tight(2, 2, 1,  [.08])

idx = [1];
plot(T_cgrid(:, idx), Y_cgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_cgrid(:, idx), X_cgrid(:, idx), 'r-.', eval_grid, out_cgrid.Z(:, idx), 'b-', eval_grid, out_cgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

[hleg1, hobj1] = legend('Sample', 'Truth', 'BHM', 'BHM 95% CI', 'Location', leglocation);

idx = [6];
plot(T_cgrid(:, idx), Y_cgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_cgrid(:, idx), X_cgrid(:, idx), 'r-.', eval_grid, out_cgrid.Z(:, idx), 'b-', eval_grid, out_cgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(a) Functional data (common grids).', 'fontsize', 16);
set(gca,'fontsize', 16)
%set(hleg1,'FontSize',12, 'color','none')
%legend boxoff


subplot_tight(2, 2, 2,  [.08])

idx = [1];
plot(T_ucgrid(:, idx), Y_ucgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_ucgrid(:, idx), X_ucgrid(:, idx), 'r-.', eval_grid, out_ucgrid.Z(:, idx), 'b-', eval_grid, out_ucgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

idx = [6];
plot(T_ucgrid(:, idx), Y_ucgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_ucgrid(:, idx), X_ucgrid(:, idx), 'r-.', eval_grid, out_ucgrid.Z(:, idx), 'b-', eval_grid, out_ucgrid.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(b) Functional data (uncommon grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 3,  [.08])
plot(eval_grid, mu_sample_cgrid, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_cgrid.mu_pgrid, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_cgrid.mu,  'b-', eval_grid, out_cgrid.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(c) Means (common grids).', 'fontsize', 16);
xlim(xlims)
ylim([-5, 5])
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.08])
plot(eval_grid, mu_sample_ucgrid, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_ucgrid.mu_pgrid, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_ucgrid.mu,  'b-', eval_grid, out_ucgrid.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(d) Means (uncommon grids).', 'fontsize', 16);
xlim(xlims)
ylim([-5, 5])
set(gca,'fontsize', 16)


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/stat_sig'))

%% plot covariance 

zlims = [-1, 8];

h = figure();

subplot_tight(2, 2, 1,  [.08,.08])
mesh(xgrid, ygrid, out_cgrid.iK)
zlim(zlims)
title('(a) BHM covariance estimate (common grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 2,  [.08,.08])
mesh(xgrid, ygrid, out_ucgrid.iK)
zlim(zlims)
title('(b) BHM covariance estimate (uncommon grids).', 'fontsize', 16);
set(gca,'fontsize', 16)


subplot_tight(2, 2, 3,  [.08,.08])
mesh(xgrid, ygrid, cov(Yfull_cgrid'))
zlim(zlims)
title('(c) Sample covariance estimate (common grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.08,.08])
mesh(xgrid, ygrid,  GausFD_ucgrid.C)
zlim(zlims)
title('(d) True covariance (common/uncommon grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/stat_cov'))


%% set up data for Nonstationary-cgrid-ucgrid plot

Y_cgrid_ns = reshape(cell2mat(GausFD_cgrid_ns.Y), length(GausFD_cgrid_ns.Y{1}), n);
X_cgrid_ns = reshape(cell2mat(GausFD_cgrid_ns.X), length(GausFD_cgrid_ns.X{1}), n);
T_cgrid_ns = reshape(cell2mat(GausFD_cgrid_ns.T), length(GausFD_cgrid_ns.T{1}), n);

Yfull_cgrid_ns = Y_cgrid_ns;
mu_sample_cgrid_ns = nanmean(Yfull_cgrid_ns, 2);


Y_ucgrid_ns = reshape(cell2mat(GausFD_ucgrid_ns.Y), length(GausFD_ucgrid_ns.Y{1}), n);
X_ucgrid_ns = reshape(cell2mat(GausFD_ucgrid_ns.X), length(GausFD_ucgrid_ns.X{1}), n);
T_ucgrid_ns = reshape(cell2mat(GausFD_ucgrid_ns.T), length(GausFD_ucgrid_ns.T{1}), n);

Yfull_ucgrid_ns = NaN(p, n); % n by p data matrix with nan's for unobserved data
    for i = 1:n
       Idx = find(ismember(pgrid, GausFD_ucgrid_ns.T{i}));
       Yfull_ucgrid_ns(Idx, i) = GausFD_ucgrid_ns.Y{i};
    end
mu_sample_ucgrid_ns = nanmean(Yfull_ucgrid_ns, 2);

%% plot stationary functional data with common and uncommon grids

leglocation='southwest';
ylims = [-15, 18]; 

h = figure();

subplot_tight(2, 2, 1,  [.08])

idx = [1];
plot(T_cgrid_ns(:, idx), Y_cgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_cgrid_ns(:, idx), X_cgrid_ns(:, idx), 'r-.', eval_grid, out_cgrid_ns.Z(:, idx), 'b-', ...
    eval_grid, out_cgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

[hleg1, hobj1] = legend('Sample', 'Truth', 'BHM', 'BHM 95% CI', 'Location', leglocation);

idx = [6];
plot(T_cgrid_ns(:, idx), Y_cgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_cgrid_ns(:, idx), X_cgrid_ns(:, idx), 'r-.', eval_grid, out_cgrid_ns.Z(:, idx), 'b-', ...
    eval_grid, out_cgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_cgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(a) Functional data (common grids).', 'fontsize', 16);
set(gca,'fontsize', 16)
%set(hleg1,'FontSize',12, 'color','none')
%legend boxoff


subplot_tight(2, 2, 2,  [.08])

idx = [1];
plot(T_ucgrid_ns(:, idx), Y_ucgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_ucgrid_ns(:, idx), X_ucgrid_ns(:, idx), 'r-.', eval_grid, out_ucgrid_ns.Z(:, idx), 'b-',...
    eval_grid, out_ucgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

idx = [10];
plot(T_ucgrid_ns(:, idx), Y_ucgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_ucgrid_ns(:, idx), X_ucgrid_ns(:, idx), 'r-.', eval_grid, out_ucgrid_ns.Z(:, idx), 'b-', ...
    eval_grid, out_ucgrid_ns.Z_CL(:, idx), 'b.', ...
    eval_grid, out_ucgrid_ns.Z_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim([-11, 11])
title('(b) Functional data (uncommon grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 3,  [.08])
plot(eval_grid, mu_sample_cgrid_ns, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_cgrid_ns.mu_pgrid, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_cgrid_ns.mu,  'b-', eval_grid, out_cgrid_ns.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(c) Means (common grids).', 'fontsize', 16);
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.08])
plot(eval_grid, mu_sample_ucgrid_ns, 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(eval_grid, GausFD_ucgrid_ns.mu_pgrid, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_ucgrid_ns.mu,  'b-', eval_grid, out_ucgrid_ns.mu_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(d) Means (uncommon grids).', 'fontsize', 16);
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/non_stat_sig'))

%%
zlims = [-1, 25];

h = figure();

subplot_tight(2, 2, 1,  [.08,.08])
mesh(xgrid, ygrid, out_cgrid_ns.iK)
zlim(zlims)
title('(a) BHM covariance estimate (common grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 2,  [.08,.08])
mesh(xgrid, ygrid, out_ucgrid_ns.iK)
zlim(zlims)
title('(b) BHM covariance estimate (uncommon grids).', 'fontsize', 16);
set(gca,'fontsize', 16)


subplot_tight(2, 2, 3,  [.08,.08])
mesh(xgrid, ygrid, cov(Yfull_cgrid_ns'))
zlim(zlims)
title('(c) Sample covariance estimate (common grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.08,.08])
mesh(xgrid, ygrid,  GausFD_ucgrid_ns.C)
zlim(zlims)
title('(d) True covariance (common/uncommon grids).', 'fontsize', 16);
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/non_stat_cov'))


%% Setup data for random grids case

Y_rgrid = reshape(cell2mat(GausFD_rgrid.Y), length(GausFD_rgrid.Y{1}), n);
X_rgrid = reshape(cell2mat(GausFD_rgrid.X), length(GausFD_rgrid.X{1}), n);
T_rgrid = reshape(cell2mat(GausFD_rgrid.T), length(GausFD_rgrid.T{1}), n);

Y_rgrid_ns = reshape(cell2mat(GausFD_rgrid_ns.Y), length(GausFD_rgrid_ns.Y{1}), n);
X_rgrid_ns = reshape(cell2mat(GausFD_rgrid_ns.X), length(GausFD_rgrid_ns.X{1}), n);
T_rgrid_ns = reshape(cell2mat(GausFD_rgrid_ns.T), length(GausFD_rgrid_ns.T{1}), n);


%% plot stationary/nonstationary functional data with random grids

leglocation='southwest';
ylims = [-11, 11]; 

h = figure();

subplot_tight(2, 2, 1,  [.08])

idx = [1];
plot(T_rgrid(:, idx), Y_rgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_rgrid(:, idx), X_rgrid(:, idx), 'r-.', eval_grid, out_rgrid.Z_cgrid(:, idx), 'b-', ...
    eval_grid, out_rgrid.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

[hleg1, hobj1] = legend('Sample', 'Truth', 'BHM', 'BHM 95% CI', 'Location', leglocation);

idx = [6];
plot(T_rgrid(:, idx), Y_rgrid(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_rgrid(:, idx), X_rgrid(:, idx), 'r-.', eval_grid, out_rgrid.Z_cgrid(:, idx), 'b-', ...
    eval_grid, out_rgrid.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim(ylims)
title('(a) Functional data (stationary).', 'fontsize', 16);
set(gca,'fontsize', 16)
%set(hleg1,'FontSize',12, 'color','none')
%legend boxoff


subplot_tight(2, 2, 2,  [.08])

idx = [1];
plot(T_rgrid_ns(:, idx), Y_rgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(T_rgrid_ns(:, idx), X_rgrid_ns(:, idx), 'r-.', eval_grid, out_rgrid_ns.Z_cgrid(:, idx), 'b-',...
    eval_grid, out_rgrid_ns.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid_ns.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)

idx = [10];
plot(T_rgrid_ns(:, idx), Y_rgrid_ns(:, idx), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
plot(T_rgrid_ns(:, idx), X_rgrid_ns(:, idx), 'r-.', eval_grid, out_rgrid_ns.Z_cgrid(:, idx), 'b-', ...
    eval_grid, out_rgrid_ns.Z_cgrid_CL(:, idx), 'b.', ...
    eval_grid, out_rgrid_ns.Z_cgrid_UL(:, idx), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim(xlims)
ylim([-10, 13])
title('(b) Functional data (nonstationary).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 3,  [.08])
plot(GausFD_rgrid.T{1}, GausFD_rgrid.Mut{1}, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_rgrid.mu_cgrid,  'b-', eval_grid, out_rgrid.mu_cgrid_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(c) Mean estimate (stationary).', 'fontsize', 16);
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.08])
plot(GausFD_rgrid_ns.T{1}, GausFD_rgrid_ns.Mut{1}, 'r-.', 'LineWidth', 2);
hold on
plot(eval_grid, out_rgrid_ns.mu_cgrid,  'b-', eval_grid, out_rgrid_ns.mu_cgrid_CI, 'b.', 'LineWidth', 2,  'MarkerSize', 10)

hold off
title('(d) Mean estimate (nonstationary).', 'fontsize', 16);
xlim(xlims)
ylim([-8, 8])
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/rand_sig'))

%%
zlims = [-1, 25];

h = figure();

subplot_tight(2, 2, 1,  [.08,.08])
mesh(xgrid, ygrid, out_rgrid.iK_cgrid)
zlim([-1, 6])
title('(a) BHM covariance estimate (stationary).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 2,  [.08,.08])
mesh(xgrid, ygrid, out_rgrid_ns.iK_cgrid)
zlim(zlims)
title('(b) BHM covariance estimate (nonstationary).', 'fontsize', 16);
set(gca,'fontsize', 16)


subplot_tight(2, 2, 3,  [.08,.08])
mesh(xgrid, ygrid, GausFD_cgrid.C)
zlim([-1, 6])
title('(c) True covariance (stationary).', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.08,.08])
mesh(xgrid, ygrid,  GausFD_cgrid_ns.C)
zlim(zlims)
title('(d) True covariance (nonstationary).', 'fontsize', 16);
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', cat(2, getenv('HOME'), '/Dropbox/FDA_Bayesian/JSS_manuscript/Figures/rand_cov'))





