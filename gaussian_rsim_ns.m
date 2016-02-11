function [simdata] = gaussian_rsim_ns(n, m, au, bu, sf, snr, nu, rho)
%Inputs
% n: # of observation, m: # of grids points for each observation
% (au, bu): domain of t ~ Uniform(au, bu)
% sf: standard deviation of the signal
% snr: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function

%Outputs
% X0 : True Signals
% 'Y': cell of signal+noise data
% 'T': cell of random grids 
% 'Mut': cell of true signal means

Ttemp = random('unif', au, bu, m, n); % random grids on t
sn = sf/ snr; % noise standard deviation

%% Stationary Gaussian data on common grid
% cells of  data

T = cell(1, n); % Cell of grids
Y = cell(1, n); % Cell of raw data
Mut = cell(1, n); % True signal mean on pooled grid
X0 = NaN(m, n); % True signal data
C = zeros(m); % covariance matrix

h = @(t)  t+1/2;

for k = 1:n

   t = sort(Ttemp(:, k));
   T{k} = t';
   t2 = t .^ 1.5 ;

   Mut{k} =  h(t) .* 3 .* sin(t2 .* 4);

   for i = 1 : m
        for j = 1 : m
            C(i, j) = h(t(i)) * h(t(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, sf^2); 
        end
   end
     
   X0(: , k) = mychol(C) * normrnd(0, 1, m, 1) + Mut{k}; % True signal, size v x n
   Y{k}  = (X0(:, k) + sn * normrnd(0, 1, m, 1))'; % Signal + noise

end

%%
simdata = struct('X0', X0, 'Y', {Y}, 'T', {T}, 'Mut', {Mut});

end




