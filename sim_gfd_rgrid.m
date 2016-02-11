%======Inputs
% n: # of observation, m: # of grids points for each observation
% sf: standard deviation of the signal
% snr: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function
% stat: generate stationary data from matern function 1 or nonstationary 0
% t ~ Uniform(au, bu)
%=======Outputs
% tau
% X0
% Y0
% 'Y', cell of noise functional data
% 'T', cell of observation grids 
% 'Xtau', true signal on the workin grid tau
% 'Mut', true mean function on pooled grid t

function [ simdata ] = sim_gfd_rgrid(n, m, au, bu, sf, snr, nu, rho, stat)

Ttemp = random('unif', au, bu, m, n); % generate random grids
sn = sf/ snr;

%% cells of  data
T = cell(1, n); % Cell of grids
Y = cell(1, n); % Cell of raw data
X = cell(1, n); % Cell of true data
Mut = cell(1, n); % True signal mean on pooled grid

C = zeros(m);
J = ones(m, 1);
h = @(t)  t+1/2;

for k = 1:n

   t = sort(Ttemp(:, k));
   T{k} = t';
  
   if stat
       Mut{k} =  3 * sin(t * 4);
       D = abs(J * t' - t * J');
       C = Matern(D, rho, nu, sf^2);
   else
       t2 = t .^ 1.5 ;
       Mut{k} =  h(t) .* 3 .* sin(t2 .* 4);
       for i = 1 : m
            for j = 1 : m
                C(i, j) = h(t(i)) * h(t(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, sf^2); 
            end
        end
   end
     
   X{k} = mychol(C) * normrnd(0, 1, m, 1) + Mut{k}; % True signal, size v x n
   Y{k}  = (X{k} + sn * normrnd(0, 1, m, 1))';

end

%%
simdata = struct('X', {X}, 'Y', {Y}, 'T', {T}, 'Mut', {Mut});

end




