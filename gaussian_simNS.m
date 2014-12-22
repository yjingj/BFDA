function [simdata] = gaussian_simNS(t, n, sf, snr, nu, rho, dense)
%Inputs
% t: pooled grid
% n: number of functional data
% sf: standard deviation of the signal
% snr: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function
% cgrid: common grid 2, or uncommon grid 0, 
% dense: 0~1, the percentage of data density

%Outputs
% 'X', cell of True signals, 
% 'Y', cell of noise functional data
% 'T', cell of grids 
% 'Idx', cell of indexes 
% 'P', number of grid points 
% 'C', ture covariance matrix on pooled grid t 
% 'mut', true mean function on pooled grid t
% 't', pooled grid


t2 = t .^ 1.5 ;
h = @(t)  t+1/2;
p = length(t);
sn = sf/ snr;

C = zeros(p); 
D = zeros(p); % Distance matrix for the pooled grid
for i = 1 : p
    for j = 1 : p
    D(i, j) = abs(t(i) - t(j));
    C(i, j) = h(t(i)) * h(t(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, sf^2); %degree of Matern Covariance Function is nu = 2.5  
    end
end



%% Stationary Gaussian data on common grid

mut = h(t) .* 3 .* sin(4 .* t2); % True signal mean on pooled grid

% cells of sparse data
P = zeros(1, n); % # of grid points for each signal
Idx = cell(1, n); % Cell of idxes of observations on the pooled grid
T = cell(1, n); % Cell of grids
X = cell(1, n); % Cell of true data
Y = cell(1, n); % Cell of raw data

% cells of common grid data
P_cgrid = zeros(1, n); % # of grid points for each signal
Idx_cgrid = cell(1, n); % Cell of idxes of observations on the pooled grid
T_cgrid = cell(1, n); % Cell of grids
X_cgrid = cell(1, n); % Cell of true data
Y_cgrid = cell(1, n); % Cell of raw data

try chol(C); catch C = topdm(C); end
%Generate functional data on pooled grid
X0 = chol(C, 'lower') * normrnd(0, 1, p, n) + repmat(mut', 1, n); % True signal, size p x n
Y0 = X0 + sn * normrnd(0, 1, p, n); % Raw signall = Truth + Noisy

%%

    % Create functional data with the common pooled grid
    for i = 1 : n
        Idx_cgrid{i} = 1:p;%sort(randsample(1:p, 0.8 * p));
        T_cgrid{i} = t(Idx_cgrid{i}); 
        P_cgrid(i) = length(Idx_cgrid{i});

        X_cgrid{i} = X0(Idx_cgrid{i}, i)'; %True signal
        Y_cgrid{i} = Y0(Idx_cgrid{i}, i)';    
    end 


    % Create functional data with uncommon grid
    for i = 1 : n
        Idx{i} = sort(randsample(1:p, ceil(dense * p))); 
        % 80% grid points are randomly sampled from the pooled grid
        T{i} = t(Idx{i}); 
        P(i) = length(Idx{i});

        X{i} = X0(Idx{i}, i)'; %True signal
        Y{i} = Y0(Idx{i}, i)';    
    end 


%%
simdata = struct('t', t, 'X0', X0, 'Y0', Y0, 'X', {X}, 'Y', {Y}, 'T', {T}, ...
    'Idx', {Idx}, 'P', P, 'X_cgrid', {X_cgrid}, 'Y_cgrid', {Y_cgrid}, 'T_cgrid', {T_cgrid}, ...
    'Idx_cgrid', {Idx_cgrid}, 'P_cgrid', P_cgrid, 'C', C, 'mut', mut, 'D', D);

end




