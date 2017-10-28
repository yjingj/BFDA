function [simdata] = sim_gfd(pgrid, n, s, r, nu, rho, dense, cgrid, stat)
% Generate functional data (common-grid, or uncommon-grid) 
% from Gaussion process with stationary Matern covariance 
% or nonstationary covariance (nonlinearly transformed Matern)

%Inputs
% pgrid: (1 X n) pooled grid
% n: number of functional curves
% s: standard deviation of the signal
% r: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function 
% dense: 0~1, the proportion of observed data on t
% cgrid: common grid 1, or uncommon grid 0,
% stat: stationary data 1, or non-stationary data 0

%Outputs
% 'Xraw_cell', cell of noisy functional data
% 'Xtrue_cell', cell of true data, 
% 'Tcell', cell of observation grids,
% 'Cov_true', ture covariance matrix on pgrid 
% 'Mean_true', true mean on pgrid

p = length(pgrid);
sn = s/ r;

if stat
    J = ones(p, 1);
    D = abs(J * pgrid - pgrid' * J'); % Distance matrix for pgrid
    Cov_true = Matern(D, rho, nu, s^2); % Matern covariance function on pgrid
    Mean_true = 3 * sin(pgrid * 4); % True signal mean on pgrid

else
    t2 = pgrid .^ 1.5 ; % transform pooled grid
    h = @(t)  t+1/2; % define transform function

    Cov_true = zeros(p); 
    for i = 1 : p
        for j = 1 : p
            Cov_true(i, j) = h(pgrid(i)) * h(pgrid(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, s^2); 
        end
    end
    Mean_true = h(pgrid) .* 3 .* sin(4 .* t2); % True signal mean on pooled grid
end

%% Gaussian functional data on common grid

Xmat_true = mychol(Cov_true) * normrnd(0, 1, p, n) + repmat(Mean_true', 1, n); % True signal, size p x n
Xmat_raw = Xmat_true + sn * normrnd(0, 1, p, n); % Raw signall = Truth + Noisy

%% 

Xraw_cell = cell(1, n); % Cell of raw data
Tcell = cell(1, n); % Cell of grids
Xtrue_cell = cell(1, n); % Cell of true data on corresponding grids

if cgrid 
    % Create functional data with the common pooled grid
    for i = 1 : n
        Tcell{i} = pgrid;
        Xtrue_cell{i} = Xmat_true(:, i)';
        Xraw_cell{i} = Xmat_raw(:, i)';         
    end 
else
    % Create functional data with uncommon grid 
    for i = 1 : n
        Idx = sort(randsample(1:p, ceil(dense * p))); 
        % partial grid points are randomly sampled from the pooled grid
        Tcell{i} = pgrid(Idx); 
        Xtrue_cell{i} = Xmat_true(Idx, i)'; % True signal on uncommon grid
        Xraw_cell{i} = Xmat_raw(Idx, i)'; % Raw signal on uncommon grid   
    end 
end

simdata = struct('Xtrue_cell', {Xtrue_cell}, 'Xraw_cell', {Xraw_cell}, 'Tcell', {Tcell}, ...
                        'Cov_true', Cov_true, 'Mean_true', Mean_true);

end




