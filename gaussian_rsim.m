function [simdata] = gaussian_rsim(n, m, au, bu, sf, snr, nu, rho)
%Inputs
% n: # of observation, m: # of grids points for each observation
% sf: standard deviation of the signal
% snr: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function
% t ~ Uniform(au, bu)

%Outputs
% tau
% X0
% Y0
% 'Y', cell of noise functional data
% 'T', cell of observation grids 
% 'Xtau', true signal on the workin grid tau
% 'Mut', true mean function on pooled grid t

Ttemp = random('unif', au, bu, m, n); % generate random grids
sn = sf/ snr;

%% Stationary Gaussian data on common grid
% cells of  data

T = cell(1, n); % Cell of grids
Y = cell(1, n); % Cell of raw data
Mut = cell(1, n); % True signal mean on pooled grid
X0 = NaN(m, n);

for k = 1:n

   t = sort(Ttemp(:, k));
   T{k} = t';
   J = ones(m, 1);
   Mut{k} =  3 * sin(t * 4);
   
   D = abs(J * t' - t * J');
   C = Matern(D, rho, nu, sf^2);
     
   X0(: , k) = mychol(C) * normrnd(0, 1, m, 1) + Mut{k}; % True signal, size v x n
   Y{k}  = (X0(:, k) + sn * normrnd(0, 1, m, 1))';

end

%%
simdata = struct('X0', X0, 'Y', {Y}, 'T', {T}, 'Mut', {Mut});

end




