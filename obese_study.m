%% load the sleep metabolic data
load('sleep.mat')
addpath(genpath(cat(2, pwd, '/mcmcdiag')))


%% Data preparation

t = 121:185;
obese_data = sleep_ee(sleep_obese == 1, t);
%non_obese_data = sleep_ee(sleep_obese == 0, t);

T = mat2cell(repmat(t, 44, 1), ones(1,44), length(t))';
Y = mat2cell(obese_data, ones(1, 44), length(t))';

%Tnon = mat2cell(repmat(t, n2, 1), ones(1,n2), length(t))';
%Ynon = mat2cell(non_obese_data, ones(1, n2), length(t))';

%BHM_non_obese = bhmcmc_fix_rho(Ynon, Tnon, 2.5, 5, 2, 1/10000, 5, 2000, 10000, 2); % for common grid
BHM_obese = bhmcmc_fix_rho(Y, T, 2.5, 5, 2, 1, 5, 2000, 10000, 2); %

save('temp.mat', 'BHM_obese')

%%
 load('obese_sp_kern.mat')
 
 %%
 h = figure(1);
subplot_tight(1, 3, 1, [0.05])

t = BHM_obese.t(1:60);
curv_oid = [15  31   40]; %for obese data
%curv_noid = [ 5 6 42]; % for non-obese data

    plot(1:length(t), BHM_obese.Yfull(1:length(t),  curv_oid),  'LineWidth', 1.5, 'Color', [0.75, 0.75, 0.75])
   hold on
   plot(1:length(t), BHM_obese.Z(1:length(t), curv_oid), 'k-', 'LineWidth', 2)
   hold on
   plot(1:length(t), BHM_obese.Z_CL(1:length(t), curv_oid), 'k-.', ...
       1:length(t), BHM_obese.Z_UL(1:length(t), curv_oid), 'k-.', 'LineWidth', 1.1)  
hold off
ylim([0.3, 1.8])
xlim([1, length(t)])
title('(a)', 'fontsize', 16);
set(gca,'fontsize', 16);


subplot_tight(1, 3, 2, [0.05])
plot(1:length(t), obese_sp(curv_oid, t)', 'k-',  'LineWidth', 2)
ylim([0.3, 1.8])
xlim([1, length(t)])
title('(b)', 'fontsize', 16);
set(gca,'fontsize', 16);


subplot_tight(1, 3, 3, [0.05])
plot(1:length(t), obese_kern(curv_oid, t)',  'k-', 'LineWidth', 2)
ylim([0.3, 1.8])
xlim([1, length(t)])
title('(c)', 'fontsize', 16);
set(gca,'fontsize', 16);



set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'obese_curv')

%% smooth data of every 50 minutes

 BHM = cell(2, 8);

for i = 1 : 8
k = (i-1) * 50 + 1;
train_id = k: (k +49);
T = mat2cell(repmat(train_id, n1, 1), ones(1,n1), length(train_id))';
Y = mat2cell(obese_data(:, train_id), ones(1, n1), length(train_id))';


Tnon = mat2cell(repmat(train_id, n2, 1), ones(1,n2), length(train_id))';
Ynon = mat2cell(non_obese_data(:, train_id), ones(1, n2), length(train_id))';

BHM_obese = bhmcmc_fix_rho(Y, T, 2.5, 5, 2, 1/100, 50, 2000, 10000, 2); %
BHM_non_obese = bhmcmc_fix_rho(Ynon, Tnon, 2.5, 5, 2, 1/100, 50, 2000, 10000, 2); % for common grid

BHM{1, i} = BHM_obese;
BHM{2, i} = BHM_non_obese;

end 
save('temp.mat', 'BHM')

%%
%curv_id = randsample(1:44, 3)



%%
T_sparse = cell(1, n1);
Y_sparse = cell(1, n1);
Tnon_sparse = cell(1, n2);
Ynon_sparse = cell(1, n2);

dense = 0.6;
for i = 1 : n1
    T_sparse{i} = sort(randsample(train_id, floor(dense * length(train_id)))); 
    Y_sparse{i} = obese_data(i, T_sparse{i});
end

for i = 1 : n2
    Tnon_sparse{i} = sort(randsample(train_id, floor(dense * length(train_id))));   
    Ynon_sparse{i} = non_obese_data(i, Tnon_sparse{i});
end

M = 10000; Burnin = 2000; nu = 2.5; delta = 5; b = 1/100; bs = 50; mat = 2; cgrid = 2;



%%

% nu = 2.5; Burnin = 2000; M = 10000;  delta = 5;  mat = 2; b = 1/100; bs = 50;
% cgrid = 2  for common grid,
BHM_obese = bhmcmc_fix_rho(Y, T, 2.5, 5, 2, 1/100, 50, 2000, 10000, 2); %
BHM_non_obese = bhmcmc_fix_rho(Ynon, Tnon, 2.5, 5, 2, 1/100, 50, 2000, 10000, 2); % for common grid

%BHM_non_obese_sparse = bhmcmc_fix_rho(Ynon_sparse, Tnon_sparse, 2.5, 30, 0, 1/100, 50, 2000, 10000, 2); 
%BHM_obese_sparse = bhmcmc_fix_rho(Y_sparse, T_sparse, 2.5, 30, 0, 1/100, 50, 2000, 10000, 2); 

cgrid = 2;
addpath(genpath(cat(2, pwd, '/PACErelease2.11')))
param_X = setOptions('regular', cgrid,'FVE_threshold',1,'corrPlot',0,'rho',-1); 
PCA_obese= FPCA(Y, T, param_X);   %perform PCA on x, common grid case
PCA_non_obese= FPCA(Ynon, Tnon, param_X);   %perform PCA on x, common grid case


% cgrid = 0, for un-common grid
%param_X = setOptions('regular', 0,'FVE_threshold',1,'corrPlot',0,'rho',-1); 
%PCA_obese_sparse= FPCA(Y, T, param_X);   %perform PCA on x, common grid case
%PCA_non_obese_sparse= FPCA(Ynon, Tnon, param_X);   %perform PCA on x, common grid case


save('obese_bhm_June22.mat', 'obese_data', 'non_obese_data', 'train_id', 'test_id',...
    'BHM_obese', 'BHM_non_obese', 'PCA_obese', 'PCA_non_obese');


!pwd | mail -s finished yjingj@GMAIL.COM

%%
save('temp.mat', 'BHM_obese', 'BHM_non_obese', 'BHM_obese_sparse', 'BHM_non_obese_sparse');

Phihat=getVal(PCA_obese,'phi');  % The estimated Eigenfunctions.
snhat2=getVal(PCA_obese,'sigma'); % The estimated noise for X(t)
lamhat=getVal(PCA_obese,'lambda'); % The estimated eigenvalues.
pcx = getVal(PCA_obese,'xi_est');  % The estimated PC scores. 

PCA_mu = getVal(PCA_obese, 'mu');
PCA_Z = pcx * Phihat';
Sigma_pace_cgrid=Phihat*diag(lamhat)*Phihat'; % The estimated covariance

