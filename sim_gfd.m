function [simdata] = sim_gfd(pgrid, n, sf, snr, nu, rho, dense, cgrid, stat)
% Generate functional data (common-grid, or uncommon-grid) 
% from Gaussion process with stationary matern covariance 
% or nonstationary covariance (nonlinearly transformed matern)

%Inputs
% pgrid: (1 X n) pooled grid
% n: number of functional curves
% sf: standard deviation of the signal
% snr: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function 
% dense: 0~1, the proportion of observed data on t
% cgrid: common grid 1, or uncommon grid 0,
% stat: stationary data 1, or non-stationary data 0

%Outputs
% 'X0', matrix of true data
% 'Y', cell of noise functional data
% 'X', cell of True data for uncommon grids, 
% 'C', ture covariance matrix on pooled grid pgrid 
% 'mu_pgrid', true mean on pooled grid pgrid

p = length(pgrid);
sn = sf/ snr;

if stat
    J = ones(p, 1);
    D = abs(J * pgrid - pgrid' * J'); % Distance matrix for pgrid
    C = Matern(D, rho, nu, sf^2); % Matern covariance function on pgrid
    mu_pgrid = 3 * sin(pgrid * 4); % True signal mean on pgrid

else
    t2 = pgrid .^ 1.5 ; % transform pooled grid
    h = @(t)  t+1/2; % define transform function

    C = zeros(p); 
    for i = 1 : p
        for j = 1 : p
            C(i, j) = h(pgrid(i)) * h(pgrid(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, sf^2); 
        end
    end
    mu_pgrid = h(pgrid) .* 3 .* sin(4 .* t2); % True signal mean on pooled grid
end

%% Gaussian functional data on common grid

X0 = mychol(C) * normrnd(0, 1, p, n) + repmat(mu_pgrid', 1, n); % True signal, size p x n
Y0 = X0 + sn * normrnd(0, 1, p, n); % Raw signall = Truth + Noisy

%% 

Y = cell(1, n); % Cell of raw data
T = cell(1, n); % Cell of grids
X = cell(1, n); % Cell of true data on corresponding grids

if cgrid 
    % Create functional data with the common pooled grid
    for i = 1 : n
        T{i} = pgrid;
        X{i} = X0(:, i)';
        Y{i} = Y0(:, i)';         
    end 
else
    % Create functional data with uncommon grid 
    for i = 1 : n
        Idx = sort(randsample(1:p, ceil(dense * p))); 
        % partial grid points are randomly sampled from the pooled grid
        T{i} = pgrid(Idx); 
        X{i} = X0(Idx, i)'; % True signal on uncommon grid
        Y{i} = Y0(Idx, i)'; % Raw signal on uncommon grid   
    end 
end

simdata = struct('X0', X0, 'X', {X}, 'Y', {Y}, 'T', {T}, ...
                        'C', C, 'mu_pgrid', mu_pgrid);

end




