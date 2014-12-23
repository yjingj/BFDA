 load('emWL.mat')
load('specZ.mat')
%% Plot raw data and smoothed data
clims = [0, 1];
curv_oid = [212    94   136];
Y = Yspecmat .^ 10;
Z = spec_BHM.Z .^ 10;
Z_sp = spec_sp .^10;
Z_kern = spec_kern .^10;

Z_UL = spec_BHM.Z_UL .^ 10;
Z_CL = spec_BHM.Z_CL .^ 10;

 %%
h = figure(1);
%curv_id = randsample(1:450, 3);
%display(curv_id);
subplot_tight(2, 2, 1,  [.05])
plot(x, Y(curv_oid, :)', 'LineWidth', 2);
title('(a)',  'fontsize', 16);
ylim([0, 16e5])
xlim([410, 700])
set(gca,'fontsize', 16,  'xticklabel',{})


subplot_tight(2, 2, 2,  [.05])
plot(xpred, Z_kern(:, curv_oid), 'LineWidth', 2);
title('(b)',  'fontsize', 16);
ylim([0, 16e5])
xlim([410, 700])
set(gca,'fontsize', 16,  'xticklabel',{})

subplot_tight(2, 2, 3,  [.05])
plot(xpred, Z_sp(:, curv_oid), 'LineWidth', 2);
title('(c)',  'fontsize', 16);
ylim([0, 16e5])
xlim([410, 700])
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.05])
plot(x, Z_CL(:, curv_oid), 'k-.', x, Z_UL(:, curv_oid), 'k-.',   'LineWidth', 1.1)
hold on
plot(x, Z(:, curv_oid), '-',  'LineWidth', 2);
hold off
title('(d)',  'fontsize', 16);
ylim([0, 16e5])
xlim([410, 700])
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1.02 1]);

print(h, '-dpdf', 'spec-curve')


%%
h = figure()
xind = 10:70;
subplot_tight(2, 2, 1,  [.05])
imagesc(x(xind), x(xind), COR(cov(Yspecmat(:, xind))), clims)
title('(a)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{})

subplot_tight(2, 2, 2,  [.05])
imagesc(x(xind), x(xind),  COR(spec_BHM.iK(xind, xind)), clims)
title('(b)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{})

subplot_tight(2, 2, 3,  [.05])
imagesc(x(xind), x(xind), COR(iK_spec_pace(xind, xind)), clims)
title('(c)', 'fontsize', 16);
set(gca,'fontsize', 16)

subplot_tight(2, 2, 4,  [.05])
imagesc(x(xind), x(xind), COR(cov(spec_BHM.Z(xind, :)')), clims)
title('(d)', 'fontsize', 16);
set(gca,'fontsize', 16)


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'spec-cor')

%% Plots for obese data
load('sleep.mat')
obese_data = sleep_ee(sleep_obese == 1, :);
non_obese_data = sleep_ee(sleep_obese == 0, :);

load('Obese_June30.mat')
 load('obese_sp_kern.mat')
 t = 81:200;

%% 
h = figure(1);
subplot_tight(3, 2, 1, [0.05])

curv_oid = [15 31   40]; %for obese data
curv_noid = [ 5 6 42]; % for non-obese data

for i = 1:8  
    BHM_obese = BHM{1, i}; 
  %  plot(BHM_obese.t, BHM_obese.mu, 'LineWidth', 2)
    %hold on
    plot(t, BHM_obese.Yfull(t,  curv_oid),  'LineWidth', 1.5, 'Color', [0.75, 0.75, 0.75])
   hold on
   plot(t, BHM_obese.Z(t, curv_oid), 'k-', 'LineWidth', 2)
   hold on
   plot(t, BHM_obese.Z_CL(t, curv_oid), 'k-.', ...
       t, BHM_obese.Z_UL(t, curv_oid), 'k-.', 'LineWidth', 1.1)  
end
hold off
ylim([0.3, 1.8])
title('(a)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{});


subplot_tight(3, 2, 2, [0.05])
for i = 1:8  
    BHM_obese = BHM{2, i}; 
  %  plot(BHM_obese.t, BHM_obese.mu, 'LineWidth', 2)
    %hold on
    plot(BHM_obese.t, BHM_obese.Yfull(:,  curv_noid),  '-.', 'LineWidth', 1.1)
   hold on
   plot(BHM_obese.t, BHM_obese.Z(:, curv_noid), '-', 'LineWidth', 2)
   hold on
   plot(BHM_obese.t, BHM_obese.Z_CL(:, curv_noid), '.', ...
       BHM_obese.t, BHM_obese.Z_UL(:, curv_noid), '.', 'MarkerSize', 6,  'Color', [0.75, 0.75, 0.75])  
end
hold off
ylim([0.3, 1.8])
title('(b)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{});


subplot_tight(3, 2, 3, [0.05])
plot(1:400, obese_sp(curv_oid, :)',   'LineWidth', 2)
ylim([0.3, 1.8])
title('(c)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{});

subplot_tight(3, 2, 4, [0.05])
plot(1:400, nonobese_sp(curv_noid, :)',   'LineWidth', 2)
ylim([0.3, 1.8])
title('(d)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{});

subplot_tight(3, 2, 5, [0.05])
plot(1:400, obese_kern(curv_oid, :)',   'LineWidth', 2)
ylim([0.3, 1.8])
title('(e)', 'fontsize', 16);
set(gca,'fontsize', 16);

subplot_tight(3, 2, 6, [0.05])
plot(1:400, nonobese_kern(curv_noid, :)',   'LineWidth', 2)
ylim([0.3, 1.8])
title('(f)', 'fontsize', 16);
set(gca,'fontsize', 16);


set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'obese_curv')

%% plot mean curve
h = figure(2);

subplot_tight(2, 2, 1, [0.05])
plot(1:400, mean(obese_data(:, 1:400), 1), 'b-',   'LineWidth', 2)
hold on
plot(1:400, mean(non_obese_data(:, 1:400), 1), 'r-',   'LineWidth', 2)
hold off
ylim([0.7, 1.3])
title('(a)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{});

subplot_tight(2, 2, 2, [0.05])

curv_noid = [ 5 6 42]; % for non-obese data
curv_oid = [15 31   40]; %for obese data

for i = 1:8  
    BHM_obese = BHM{1, i}; 
   plot(BHM_obese.t, BHM_obese.mu, 'b-', 'LineWidth', 2)
    hold on
    plot(BHM_obese.t, BHM_obese.mu_CI, 'b.', 'MarkerSize', 6)
end

for i = 1:8  
    BHM_nonobese = BHM{2, i}; 
   plot(BHM_nonobese.t, BHM_nonobese.mu, 'r-', 'LineWidth', 2)
   hold on
    plot(BHM_nonobese.t, BHM_nonobese.mu_CI, 'r.', 'MarkerSize', 6)
end
hold off

ylim([0.7, 1.3])
title('(b)', 'fontsize', 16);
set(gca,'fontsize', 16,  'xticklabel',{});



subplot_tight(2, 2, 3, [0.05])
plot(1:400, mean(obese_sp, 1),  'b-', 'LineWidth', 2)
hold on
plot(1:400, mean(nonobese_sp, 1), 'r-',  'LineWidth', 2)
hold off
ylim([0.7, 1.2])
title('(c)', 'fontsize', 16);
set(gca,'fontsize', 16);

subplot_tight(2, 2, 4, [0.05])
plot(1:400, mean(obese_kern, 1), 'b-',  'LineWidth', 2)
hold on
plot(1:400, mean(nonobese_kern, 1), 'r-',  'LineWidth', 2)
hold off
ylim([0.7, 1.3])
title('(d)', 'fontsize', 16);
set(gca,'fontsize', 16);

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'obese_mean')








