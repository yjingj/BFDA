%======Inputs
% n: # of observation, m: # of grid points for each observation
% s: signal standard deviation 
% r: signal to noise ratio
% nu: smoothness order in the Matern function
% rho: spatio scale parameter in the Matern function
% stat: generate data from the stationary GP with value 1 
%       or nonstationary GP with value 0
% t ~ Uniform(au, bu)

%=======Outputs
% 'Xtrue_cell', cell of true functional data
% 'Xraw_cell', cell of noise functional data
% 'Tcell', cell of observation grids 
% 'Mean_true_cell', true mean function on pooled grid t

function [ simdata ] = sim_gfd_rgrid(n, m, au, bu, s, r, nu, rho, stat)

Ttemp = random('unif', au, bu, m, n); % generate random grids
sn = s/ r;

%% cells of  data
Tcell = cell(1, n); % Cell of grids
Xraw_cell = cell(1, n); % Cell of raw data
Xtrue_cell = cell(1, n); % Cell of true data
Mean_true_cell = cell(1, n); % True signal mean on pooled grid

C = zeros(m);
J = ones(m, 1);
h = @(t)  t+1/2;

for k = 1:n

   t = sort(Ttemp(:, k));
   Tcell{k} = t';
  
   if stat
       Mean_true_cell{k} =  3 * sin(t * 4);
       D = abs(J * t' - t * J');
       C = Matern(D, rho, nu, s^2);
   else
       t2 = t .^ 1.5 ;
       Mean_true_cell{k} =  h(t) .* 3 .* sin(t2 .* 4);
       for i = 1 : m
            for j = 1 : m
                C(i, j) = h(t(i)) * h(t(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, s^2); 
            end
        end
   end
     
   Xtrue_cell{k} = mychol(C) * normrnd(0, 1, m, 1) + Mean_true_cell{k}; % True signal, size v x n
   Xraw_cell{k}  = (Xtrue_cell{k} + sn * normrnd(0, 1, m, 1))';

end

%%
simdata = struct('Xtrue_cell', {Xtrue_cell}, 'Xraw_cell', {Xraw_cell}, 'Tcell', {Tcell}, 'Mean_true_cell', {Mean_true_cell});

end




