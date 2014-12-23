close all
clear all
load('station.mat')

%% Compare with lowess
[n, p] = size(Y);
Z3 = zeros(n, p);
for i = 1 : n
  Z3(i, :) = smooth( Y(i, :), 0.4, 'lowess');
end

%Local regression using weighted linear least squares 
%and a 1st degree polynomial model

Z4 = zeros(n, p);
for i = 1 : n
  [SP  Z4(i, :)] = spaps(x, Y(i, :), 0.9);
end
%B-form of the smoothest function f that lies within
%the given tolerance tol of the given data points 

plot(Y(1, :), 'k')
hold on
plot(Z3(1, :), 'g')
plot(Z4(1, :), 'm')
plot(Z(:, 1), 'r')
hold off

% Display the MSE

display(['MSE for raw signal', num2str(mean(mean((Y - X).^2)))])
  % MSE for raw signal
display(['MSE for smoothed signal by B-form splines functions: ', num2str(mean(mean((Z4 - X).^2)))])
 %MSE for smoothed signal by B-form splines functions
display(['MSE for smoothed signal by lowess: ', num2str(mean(mean((Z3 - X).^2)))])
%MSE for smoothed signal by lowess
display(['MSE for smoothed signal by MCMC: ', num2str(mean(mean((Z' - X).^2)))])

display(['MSE for mean signal by MCMC: ', num2str(mean((mu - mut(:, 1)).^2))])
 % MSE for smoothed signal mean
 display(['MSE for mean signal by B-splines: ', num2str(mean((mean(Z4)' - mut(:, 1)).^2))])
 % MSE for mean signal by B-splines
display(['MSE for mean signal by lowess: ', num2str(mean((mean(Z3)' - mut(:, 1)).^2))])
 % MSE for mean signal by lowess


%%
% Sum square error

 % MSE for smoothed signal by MCMC
 display(['MSE for sample estimated COV: ', num2str(mean(mean((C - cov(Y)).^2)))])
 % MSE for sample estimated COV
display(['MSE for MCMC estimated COV: ', num2str(mean(mean((C - iK).^2)))])
 % MSE for MCMC estimated COV

% Sum square error

%%
close all
h = figure(4);
subplot(2, 1, 1)
plot(x, Y(1:10, :), 'LineWidth', 2);
%title(char('First 10 raw stationary functional data curves'), 'fontsize', 16);
xlabel('t','fontsize', 16)
set(gca,'fontsize', 16)

subplot(2, 1, 2)
plot(x, Z(:, 1:10), 'LineWidth', 2);
%title(char('First 10 Bayesian smoothed functional data curves'), 'fontsize', 16);
xlabel('t','fontsize', 16)
set(gca,'fontsize', 16)
set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'fdCurves-obs')


%% Plot noise standard deviation
h = figure(1);
subplot(2, 1, 1)
plot(1./sqrt(TaoOut))
title('Trace plot of noise standard deviation (First group of MCMC)','fontsize', 16)
set(gca,'fontsize', 16)
subplot(2, 1, 2)
hist(1./sqrt(TaoOut))
title('Histogram of noise standard deviation (First group of MCMC)','fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'sd-noise1')

%% Plot noise standard deviation
h = figure(1);
subplot(2, 1, 1)
plot(1./sqrt(TaoOut))
title('Trace plot of \sigma_{\epsilon}','fontsize', 16)
set(gca,'fontsize', 16)
subplot(2, 1, 2)
hist(1./sqrt(TaoOut))
title('Histogram of \sigma_{\epsilon}','fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'sd-noise2')

%% Plot mean function
musort = sort(muOut, 2);

h = figure(2);
plot(x, mut(:, 1), 'K--', x, mean(Y), 'b.', x, mu, 'r-', x, musort(:, 9750), 'r-.', ...
x, musort(:, 250), 'r-.',  'MarkerSize',20, 'LineWidth', 3)
    ylim([-2, 3.5])
    set(gca,'fontsize', 16)
%title('Mean Trajectories', 'FontSize', 16)
xlabel('t', 'FontSize', 16)
[hleg1, hobj1] = legend('True GP Mean', 'Observed Signal Mean', 'Bayesian Estimated Mean', ...
    '95% Credible Interval', 'Location', 'Best');
set(hleg1, 'position', [.4 .25 0.01 .1],  'FontSize', 25) 

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'mean-curve')


%% Plot correlation matrix
[n, p] = size(Y);
SC = cov(Y); % Sample covariance function
diag_SC = diag(SC);
SCR = SC ./ (sqrt(diag_SC) * sqrt(diag_SC)' );

diag_iK = diag(iK);
COR = iK ./ (sqrt(diag_iK) * sqrt(diag_iK)' );

diag_C = diag(C);
Cor_true = C ./  ((sqrt(diag_C) * sqrt(diag_C)' ));
clims = [ min([min(SCR), min(C ./ 5), min(COR)]),  1 ];

%%
h = figure(1);
imagesc(C./ 5, clims);
colorbar
%title('Heat Map of the True Correlation Matrix.', 'FontSize', 16)
set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
set(gca,'fontsize', 16)
print(h, '-dpdf', 'heat-true')

%%

h = figure(2);
imagesc(SCR, clims);
colorbar
%title('Heat Map of the Sample Estimated Corrlation Matrix.', 'FontSize', 16)
set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
set(gca,'fontsize', 16)
print(h, '-dpdf', 'heat-sample')
%%
h = figure(3);
imagesc(COR, clims);
colorbar
%title('Heat Map of the Bayesian Estimated Corrlation Matrix.', 'FontSize', 16)
set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
set(gca,'fontsize', 16)
print(h, '-dpdf', 'heat-mcmc')

%% Plot 4 data curves
SSE = findsse(X, Z');
[SSE1, Idx] = sort(SSE);
selectIdx = Idx([p, ceil(p/3 * 2), ceil(p/3),  1]);

%%
Zsort = sort(ZOut, 3);
h = figure(4);
for i = 1:4
    j = selectIdx(i);
    subplot(2, 2, i)
    plot(x, X(j, :), 'k--', x, Y(j, :), 'b.', x, Z(:, j), 'r', ...
        x, Zsort(:, j, 9750), 'r-.', x, Zsort(:, j, 250), 'r-.', 'MarkerSize',20,  'LineWidth', 3)
    xlabel('t', 'FontSize', 16)
    title(['Signal ', num2str(j)],  'FontSize', 16)
    set(gca,'fontsize', 16)
    if i == 4
[hleg2, hobj2]  =   legend('True Signal', 'Observed Signal', 'Bayesian Smoothed Signal',  ...
         '95% Credible Interval', 'Location', 'Best');
     set(hleg2, 'position', [.65 .15 .01 .01]) 
    end
end

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'curves')

%% Plot diagnal variances
iKsort = sort(iKOut, 3);

 h= figure(5);
plot(1:p, diag(C), 'k--', 1:p, diag(SC), 'b.', 1:p, diag(iK), 'r', ...
    1:p, diag(iKsort(:, :, 9750)), 'r-.', 1:p, diag(iKsort(:, :, 250)), 'r-.', ...
     'MarkerSize',20, 'LineWidth', 3)
%title('Variance of Each Dimension', 'fontsize', 16)

[hleg3, hobj3] = legend('True Variance', 'Sample Estimated Variance',...
    'Bayesian Estimated Variance', '95% Credible Interval', 'Location', 'Best');
set(hleg3,  'position', [.3 .75 0.01 .1], 'FontSize', 25) 

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
set(gca,'YLim', [3  8.5]);
set(gca,'fontsize', 16)
print(h, '-dpdf', 'sim-var')

%% Plot diagnal and selected rows of covariance matrix

[SSE2, Idx2] = sort(findsse(C, iK));

selectIdx2 = Idx2([p, ceil(p/3 * 2), ceil(p/3),  1]);

%% Plot diagnal and selected rows of covariance matrix
h = figure(9);
subplot(2, 2, 1)
for i = 1 : 4
 j = selectIdx2(i);
subplot(2, 2,  i)
plot(1:p, C(j, :), 'k-.', 1:p, SC(j, :), 'b--', 1:p, iK(j, :),...
    'r', 'LineWidth', 2)
title(['Row ', num2str(j)],'fontsize', 16)
set(gca,'fontsize', 16)
end

[hleg3, hobj3]  =   legend('True Cov-Matrix', 'Sample Estimated Cov-Matrix',...
    'Bayesian Estimated Cov-Matrix', 'Location', 'Best');
set(hleg3, 'position', [.4 .33 0.01 .1]) 

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'cov-rows')


%% Trace plot of log proportional joint likelihood function 
h = figure(6);
plot(loglike)
title('Trace plot of log proportional joint-likelihood', 'fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'loglike-trace')

%%
h = figure(7);
subplot(2, 1, 1)
plot(rsOut)
title('Trace plot of \sigma_s^2', 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(2, 1, 2)
hist(rsOut)
title('Histogram of \sigma_s^2', 'fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'sigma-s2')

%%
h = figure(8);
subplot(2, 1, 1)
plot(rtOut)
title('Trace plot of \gamma_{\theta_1}', 'fontsize', 16)
set(gca,'fontsize', 16)
subplot(2, 1, 2)
hist(rtOut)
title('Histogram of \gamma_{\theta_1}', 'fontsize', 16)
set(gca,'fontsize', 16)

set(h, 'PaperOrientation','landscape');
set(h, 'PaperUnits','normalized');
set(h, 'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', 'gamma-theta1')
