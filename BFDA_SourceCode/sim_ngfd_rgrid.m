%======Inputs
% n: # of observation, m: # of grids points for each observation
% sf: standard deviation of the signal
% snr: signal to noise ratio
% nu: order in matern function
% rho: spacial scale parameter in matern function
% stat: generate stationary data from matern function 1 or nonstationary 0
% t ~ Uniform(au, bu)

%=======Outputs
% 'Xtrue_cell', cell of true functional data
% 'Xraw_cell', cell of noise functional data
% 'Tcell', cell of observation grids 
% 'Mean_true_cell', true mean function on pooled grid t

function [ simdata ] = sim_ngfd_rgrid(n, p, au, bu, sf, snr, nu, rho, stat)

Ttemp = random('unif', au, bu, p, n); % generate random grids
sn = sf/ snr;

%% cells of  data
Tcell = cell(1, n); % Cell of grids
Xraw_cell = cell(1, n); % Cell of raw data
Xtrue_cell_gp = cell(1, n); % Cell of true Gaussian data
Xtrue_cell = cell(1, n); % Cell of true Non-Gaussian data
Mean_true_cell = cell(1, n); % True signal mean on pooled grid

C = zeros(p);
J = ones(p, 1);
h = @(t)  t+1/2;



for k = 1:n

   t = sort(Ttemp(:, k));
   Tcell{k} = t';
  
   if stat
       Mean_true_cell{k} = 3 * sin(t * 4);
       D = abs(J * t' - t * J');
       C = Matern(D, rho, nu, sf^2);
   else
       t2 = t .^ 1.5 ;
       Mean_true_cell{k} =  h(t) .* 3 .* sin(t2 .* 4);
       for i = 1 : p
            for j = 1 : p
                C(i, j) = h(t(i)) * h(t(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, sf^2); 
            end
        end
   end
     
   Xtrue_cell_gp{k} = mychol(C) * normrnd(0, 1, p, 1) + Mean_true_cell{k}; % True signal, size v x n
   Xtrue_cell{k} = 0.2 .* ((Xtrue_cell_gp{k}).^(2) - 1) + Xtrue_cell_gp{k};
   Xraw_cell{k}  = (Xtrue_cell{k} + sn * normrnd(0, 1, p, 1))';

end

%%
simdata = struct('Xtrue_cell', {Xtrue_cell}, 'Xraw_cell', {Xraw_cell}, ...
    'Xtrue_cell_gp', {Xtrue_cell_gp}, ...
    'Tcell', {Tcell}, 'Mean_true_cell', {Mean_true_cell});

end




