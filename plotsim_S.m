% close all
% load('StationaryGauss.mat')

%% Stationary Gaussian Data
t = GausData.t;
C = GausData.C;
D = GausData.D;
X0 = GausData.X0;
Y0 = GausData.Y0;
mut = GausData.mut;

X = GausData.X;
Y = GausData.Y;
T = GausData.T;

X_cgrid = GausData.X_cgrid;
Y_cgrid = GausData.Y_cgrid;
T_cgrid = GausData.T_cgrid;


%%
Z = BHM_cgrid.Z;  Z_CL = BHM_cgrid.Z_CL; Z_UL = BHM_cgrid.Z_UL;
iK = BHM_cgrid.iK;
mu = BHM_cgrid.mu;
Yfull = BHM_cgrid.Yfull;

Z_sparse = BHM_sparse.Z; Zsparse_CL = BHM_sparse.Z_CL; Zsparse_UL = BHM_sparse.Z_UL;
iK_sparse = BHM_sparse.iK;
mu_sparse = BHM_sparse.mu;
Yfull_sparse = BHM_sparse.Yfull;

%% coverage probability

Covprob(mut', BHM_cgrid.mu_CI(:, 1), BHM_cgrid.mu_CI(:, 2), 1)
Covprob(X0, BHM_cgrid.Z_CL, BHM_cgrid.Z_UL, 2)
Covprob(C, BHM_cgrid.iK_CL, BHM_cgrid.iK_UL, 2)

Covprob(mut', BHM_sparse.mu_CI(:, 1), BHM_sparse.mu_CI(:, 2), 1)
Covprob(X0, BHM_sparse.Z_CL, BHM_sparse.Z_UL, 2)
Covprob(C, BHM_sparse.iK_CL, BHM_sparse.iK_UL, 2)




%%
addpath(genpath(cat(2, pwd, '/PACErelease2.11')))

Phihat=getVal(PCA_cgrid,'phi');  % The estimated Eigenfunctions.
snhat2=getVal(PCA_cgrid,'sigma'); % The estimated noise for X(t)
lamhat=getVal(PCA_cgrid,'lambda'); % The estimated eigenvalues.
pcx = getVal(PCA_cgrid,'xi_est');  % The estimated PC scores. 

PCA_mu = getVal(PCA_cgrid, 'mu');
PCA_Z = pcx * Phihat';
Sigma_pace_cgrid=Phihat*diag(lamhat)*Phihat'; % The estimated covariance


pcx_sparse = getVal(PCA_sparse, 'xi_est');  % The estimated PC scores. 
Phihat_sparse=getVal(PCA_sparse,'phi');  % The estimated Eigenfunctions.
snhat2_sparse=getVal(PCA_sparse,'sigma'); % The estimated noise for X(t)
lamhat_sparse=getVal(PCA_sparse,'lambda'); % The estimated eigenvalues.

PCA_mu_sparse = getVal(PCA_sparse, 'mu');
PCA_Z_sparse = pcx_sparse * Phihat_sparse';
Sigma_pace_sparse=Phihat_sparse*diag(lamhat_sparse)*Phihat_sparse'; % The estimated covariance



%%
rimse(t, X0, BHM_cgrid.Z, 0)
rimse(t, mut, BHM_cgrid.mu', 0)
rimse(t, C, BHM_cgrid.iK, 2)

rmse(X0, Z)
rmse(C, iK)
rmse(mut, mu')

rmse(X0, Z_sparse)
rmse(C, iK_sparse)
rmse(mut, mu_sparse')


%%
clims = [0, 1];

 h = figure();
subplot_tight(3, 2, 1,  [.05])
plot(t, Yfull(:, [1 25]), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(t, Z(:, [1 25]),  'k-', t, Z_CL(:, [1, 25]), 'k-.', t, Z_UL(:, [1, 25]), 'k-.', t, X0(:, [1, 25]), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim([0, pi/2])
ylim([-8, 10])
title('(a)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 2,  [.05])
plot(t, Yfull_sparse(:, [1, 25]), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(t, Z_sparse(:, [1 25]),  'k-', t, Zsparse_CL(:, [1, 25]), 'k-.', t, Zsparse_UL(:, [1, 25]), 'k-.', t, X0(:, [1, 25]), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
title('(b)', 'fontsize', 16);
xlim([0, pi/2])
ylim([-8, 10])
set(gca,'fontsize', 16)

subplot_tight(3, 2, 3,  [.05])
imagesc(t, t, COR(iK), clims)
title('(c)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 4,  [.05])
imagesc(t, t, COR(iK_sparse), clims)
title('(d)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 5,  [.05])
imagesc(t, t, COR(cov(Yfull')), clims)
title('(e)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 6,  [.05])
imagesc(t, t, COR(C), clims)
title('(f)', 'fontsize', 16);
set(gca,'fontsize', 16)


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', 'curves')


%% 








%%

[xgrid,ygrid]=meshgrid(GausData.t); % 3D plots of the covariance matrices

h=figure();
subplot(2, 2, 1)
mesh(xgrid, ygrid, BHM_cgrid.iK)
zlim([0, 5])
title('Covariance estimate with BHM', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(2, 2, 2)
mesh(xgrid, ygrid, Tprior_cgrid.iK)
zlim([0, 5])
title('Covariance estimate with t-folded prior', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(2, 2, 3)
mesh(xgrid, ygrid, Matern_cgrid.iK)
zlim([0, 5])
title('Covariance estimate with Matern model', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(2, 2, 4)
mesh(xgrid,ygrid,C)
zlim([0, 5])
title('True covariance', 'fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);

print(h, '-dpdf', 'matern_cov')

%%
h = figure();
subplot(3, 2, 1)
plot(t, Y0(:, 1:10), 'LineWidth', 2);
xlim([0, pi/2])
ylim([-8, 10])
title(char('Functional data observed on common grid (non-sparse)'), 'fontsize', 16);
xlabel('(a)','fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 2)
plot(t, BHM_cgrid.Z(:, 1:10), 'LineWidth', 2);
title(char('Signal estimates from BHM'), 'fontsize', 16);
xlim([0, pi/2])
ylim([-8, 10])
xlabel('(c)','fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 3)
plot(t, Tprior_cgrid.Z(:, 1:10), 'LineWidth', 2);
title(char('Signal estimates with t-folded prior'), 'fontsize', 16);
xlim([0, pi/2])
ylim([-8, 10])
xlabel('(d)','fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 4)
plot(t, Matern_cgrid.Z(:, 1:10), 'LineWidth', 2);
title(char('Signal estimates from Matern model'), 'fontsize', 16);
xlim([0, pi/2])
ylim([-8, 10])
xlabel('(d)','fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 5)
plot(t, X0(:, 1:10), 'LineWidth', 2);
title(char('True functional data'), 'fontsize', 16);
xlim([0, pi/2])
ylim([-8, 10])
xlabel('(e)','fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 6)
plot(t, mut, 'K--', t, BHM_cgrid.mu, 'b-', t, Tprior_cgrid.mu, 'r-.', ...
t, Matern_cgrid.mu, 'g-.', 'MarkerSize',20, 'LineWidth', 3)
title(char('Functional mean'), 'fontsize', 16);
xlim([0, pi/2])
ylim([-8, 10])
xlabel('(f)','fontsize', 16)
%set(gca,'fontsize', 16)
[hleg1, hobj1] = legend('True GP Mean', 'BHM',...
    't-folded prior', 'Matern model', 'position', 'Best');
set(hleg1, 'position', [.82 .22 0.01 .1],  'FontSize', 10) 
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'matern_curves')