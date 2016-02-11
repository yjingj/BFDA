function [C] = cov_ns(T, sf, nu, rho)
%Inputs
% T: a random grid
% sf: standard deviation of the signal
% nu: order in matern function
% rho: spacial scale parameter in matern function

%Outputs
%C: non-stationary covaricance

m = length(T); 
C = zeros(m); % covariance matrix

h = @(t)  t+1/2;

t = sort(T);

t2 = t .^ 1.5 ;

for i = 1 : m
   for j = 1 : m
       C(i, j) = h(t(i)) * h(t(j)) * Matern(abs(t2(i) - t2(j)), rho, nu, sf^2); 
   end
end
     
end




