function r=rimse(t, data, estimate, twoD)
% Function to calculate root integrated mean square error from a function
% estimate over its domain and the corresponding estimates.
% data has size n * p
% Note: data and estimates have to be vectors of same size, over grid t

% delete records with NaNs in both datasets first

mse = (data - estimate).^2; 

if (~twoD)
    
    r = sqrt( mean(trapz(t, mse) ) );
    
else 
    
    r = sqrt(trapz(t, trapz(t, mse)));
    
end

