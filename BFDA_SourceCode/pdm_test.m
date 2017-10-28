% pdm_test() function is used for goodness-of-fit diagnostics
% w: a nxKxJ matrix of standardized residuals, with K groups, 
%       J mcmc iterations, n residuals
% L: number of cells, default value is 4
% M: number of order statistics used for test, default is 100
function[pmin, pdm, df] = pdm_test(w, L, M)
    
    if isempty(w)
        error('Require input of residuals!');
    end
    [n, K, J] = size(w);
    
    if isempty(L)
        L = 4;
    end
    
    if isempty(M)
        M = 100;
    elseif M > J;
        M = J;
        display('M should not be larger than the third dimension of w!');
    end
  
    df = K * (L - 1); % specify degrees of freedom for chi-square dist.
    q_L = linspace(0.25, 0.75, L-1); % 
    split_L = norminv(q_L, 0, 1); % split quantiles length L-1
    p_L = [q_L, 1] - [0, q_L]; % cumulative density per cell
    e_L = n .* p_L; % expected number of observations
    
    pdm = zeros(1, J);
    
    for j = 1:J
        for k = 1:K
            w_k_j = sort(reshape(w(:, k, j), 1, n));
            for l = 1:L
                if l == 1
                    w_k_l = w_k_j(w_k_j < split_L(l));
                elseif l < L
                    w_k_l = w_k_j(w_k_j < split_L(l) & w_k_j > split_L(l-1));
                else
                    w_k_l = w_k_j(w_k_j > split_L(L-1));
                end
                o_k_l = length(w_k_l);
                pdm(j) = pdm(j) + ((o_k_l - e_L(l)) / sqrt(e_L(l)))^2;   
            end
        end
    end
    
    pdm = sort(pdm);
    q_M = linspace(0.005, 0.995, M); % specify test quantiles
    order_pdm = floor(quantile(1:J, q_M));
    q_pdm = pdm(order_pdm);
    pval = min(1, (J .* (1 - chi2cdf(q_pdm, df))) ./ (J - order_pdm + 1) );
    pmin = min(pval);

end



