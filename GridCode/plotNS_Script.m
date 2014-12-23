% close all
%load('StationaryGauss_NS.mat')

%% Stationary Gaussian Data
tns = GausData_NS.t;
Cns = GausData_NS.C;
Dns = GausData_NS.D;
X0ns = GausData_NS.X0;
Y0ns = GausData_NS.Y0;
mutns = GausData_NS.mut;


Xns = GausData_NS.X;
Yns = GausData_NS.Y;
Tns = GausData_NS.T;

Xns_cgrid = GausData_NS.X_cgrid;
Yns_cgrid = GausData_NS.Y_cgrid;
Tns_cgrid = GausData_NS.T_cgrid;

%%
Zns = BHMns_cgrid.Z; Zns_CL = BHMns_cgrid.Z_CL; Zns_UL = BHMns_cgrid.Z_UL;
iKns = BHMns_cgrid.iK;
muns = BHMns_cgrid.mu;
Yfullns = BHMns_cgrid.Yfull;

Zns_sparse = BHMns_sparse.Z; Zsparse_CLns = BHMns_cgrid.Z_CL; Zsparse_ULns = BHMns_cgrid.Z_UL;
iKns_sparse = BHMns_sparse.iK;
muns_sparse = BHMns_sparse.mu;
Yfullns_sparse = BHMns_sparse.Yfull;
%%

Covprob(X0ns, BHMns_cgrid.Z_CL, BHMns_cgrid.Z_UL, 2)
Covprob(mutns', BHMns_cgrid.mu_CI(:, 1), BHMns_cgrid.mu_CI(:, 2), 1)
Covprob(Cns, BHMns_cgrid.iK_CL, BHMns_cgrid.iK_UL, 2)


Covprob(X0ns, BHMns_sparse.Z_CL, BHMns_sparse.Z_UL, 2)
Covprob(mutns', BHMns_sparse.mu_CI(:, 1), BHMns_sparse.mu_CI(:, 2), 1)
Covprob(Cns, BHMns_sparse.iK_CL, BHMns_sparse.iK_UL, 2)

%%

addpath(genpath(cat(2, pwd, '/PACErelease2.11')))

Phihat=getVal(PCAns_cgrid,'phi');  % The estimated Eigenfunctions.
snhat2=getVal(PCAns_cgrid,'sigma'); % The estimated noise for X(t)
lamhat=getVal(PCAns_cgrid,'lambda'); % The estimated eigenvalues.
pcx = getVal(PCAns_cgrid,'xi_est');  % The estimated PC scores. 

PCAns_mu = getVal(PCAns_cgrid, 'mu');
PCAns_Z = pcx * Phihat';
Sigmans_pace_cgrid=Phihat*diag(lamhat)*Phihat'; % The estimated covariance

pcx_sparse = getVal(PCAns_sparse, 'xi_est');  % The estimated PC scores. 
Phihat_sparse=getVal(PCAns_sparse,'phi');  % The estimated Eigenfunctions.
snhat2_sparse=getVal(PCAns_sparse,'sigma'); % The estimated noise for X(t)
lamhat_sparse=getVal(PCAns_sparse,'lambda'); % The estimated eigenvalues.

PCAns_mu_sparse = getVal(PCAns_sparse, 'mu');
PCAns_Z_sparse = pcx_sparse * Phihat_sparse';
Sigmans_pace_sparse=Phihat_sparse*diag(lamhat_sparse)*Phihat_sparse'; % The estimated covariance



%%
rimse(t, X0ns, BHMns_cgrid.Z, 0)
rimse(t, mutns, BHMns_cgrid.mu', 0)
rimse(t, Cns, BHMns_cgrid.iK, 2)

rmse(X0ns, Zns)
rmse(Cns, iKns)
rmse(mutns, muns')

rmse(X0ns, Zns_sparse)
rmse(Cns, iKns_sparse)
rmse(mutns, muns_sparse')


%% 
 
clims = [0, 1];

 h = figure();
 
 curv_id = [1, 9];
 
subplot_tight(3, 2, 1,  [.05])
plot(tns, Yfullns(:, curv_id), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(tns, Zns(:, curv_id),  'k-', tns, Zns_CL(:, curv_id), 'k-.', t, Zns_UL(:, curv_id), 'k-.', t, X0ns(:, curv_id), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
xlim([0, pi/2])
ylim([-11, 16])
title('(a)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 2,  [.05])
plot(tns, Yfullns_sparse(:, curv_id), 'LineWidth', 2, 'Color', [0.75, 0.75, 0.75]);
hold on
plot(tns, Zns_sparse(:, curv_id),  'k-', tns, Zsparse_CLns(:, curv_id), 'k-.', t, Zsparse_ULns(:, curv_id), 'k-.', t, ...
    X0ns(:, curv_id), 'b.', 'LineWidth', 2,  'MarkerSize', 10)
hold off
title('(b)', 'fontsize', 16);
xlim([0, pi/2])
ylim([-11, 16])
set(gca,'fontsize', 16)

subplot_tight(3, 2, 3,  [.05])
imagesc(t, t, COR(iKns), clims)
title('(c)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 4,  [.05])
imagesc(t, t, COR(iKns_sparse), clims)
title('(d)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 5,  [.05])
imagesc(t, t, COR(cov(Yfullns')), clims)
title('(e)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(3, 2, 6,  [.05])
imagesc(t, t, COR(Cns), clims)
title('(f)', 'fontsize', 16);
set(gca,'fontsize', 16)


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', 'non-curves')



%% Plot covariance matrix

SCov = cov(Y0');  % PACE estimated covariance matrix
diag_SC = diag(SCov);
SCR = SCov ./ (sqrt(diag_SC) * sqrt(diag_SC)' );

R = C ./ (sqrt(diag(C)) * sqrt(diag(C))' );
BHM_R = BHM_cgrid.iK ./ (sqrt(diag(BHM_cgrid.iK)) * sqrt(diag(BHM_cgrid.iK))' );
diag_iK = diag(iKns_sparse); % Bayesian estimate
R_sparse = iKns_sparse ./ (sqrt(diag_iK) * sqrt(diag_iK)' );
pace_R = Sigmans_pace_cgrid ./ (sqrt(diag(Sigmans_pace_cgrid)) * sqrt(diag(Sigmans_pace_cgrid))' );

[xgrid,ygrid]=meshgrid(t); % 3D plots of the covariance matrices

figure()
mesh(xgrid,ygrid, Sigmans_pace_cgrid)
figure()
mesh(xgrid, ygrid, BHM_cgrid.iK)

%%
clims = [0, 1];

h=figure();
subplot(3, 2, 1)
imagesc(t, t, iKns ./ (sqrt(diag(iKns)) *  sqrt(diag(iKns))'), clims)
colorbar
title('Bayesian correlation estimate (non-sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 2)
imagesc(t, t, iKns_sparse ./ (sqrt(diag(iKns_sparse)) *  sqrt(diag(iKns_sparse))'), clims)
colorbar
title('Bayesian correlation estimate (sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 3)
imagesc(t, t, Sigmans_pace_cgrid ./ (sqrt(diag(Sigmans_pace_cgrid)) *  sqrt(diag(Sigmans_pace_cgrid))'), clims)
colorbar
title('PACE correlation estimate (non-sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 4)
imagesc(t, t, Sigmans_pace_sparse ./ (sqrt(diag(Sigmans_pace_sparse)) *  sqrt(diag(Sigmans_pace_sparse))'), clims)
colorbar
title('PACE correlation estimate (sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 5)
imagesc(t, t, SCov ./ (sqrt(diag(SCov)) *  sqrt(diag(SCov))'), clims) 
colorbar
title('Sample correlation (non-sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 6)
imagesc(t, t, R, clims);
colorbar
title('True correlation', 'fontsize', 16)
set(gca,'fontsize', 16)


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'NS_cor_heat')

%%
h=figure();
subplot(3, 2, 1)
mesh(xgrid, ygrid, iKns)
zlim([0, 18])
title('Bayesian covariance estimate (non-sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 2)
mesh(xgrid, ygrid, iKns_sparse)
zlim([0, 18])
title('Bayesian covariance estimate (sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 3)
mesh(xgrid, ygrid, Sigmans_pace_cgrid)
zlim([0, 18])
title('PACE covariance estimate (non-sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 4)
mesh(xgrid, ygrid, Sigmans_pace_sparse)
zlim([0, 18])
title('PACE covariance estimate (sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 5)
mesh(xgrid,ygrid, SCov)
zlim([0, 18])
title('Sample covariance estimate (non-sparse)', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(3, 2, 6)
mesh(xgrid,ygrid,C)
zlim([0, 18])
title('True covariance', 'fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);

 print(h, '-dpdf', 'NS_cov')

%%

