%This function obtain the marginal log-likelihood for AIC1 or BIC1
%It is used in no_AIC or no_BIC functions.
%
%[invalid, logLik1] = getLogLik1(y,lambda, muSub, phiSub, sigma,regular)
% Input:
%     y: 1*n cell array of observed repeated measurements for n subjects
%  lambda: 1*k vector of eigenvalues
%  muSub: 1*n cell array of mean functions evaluated at the t{i}'s
% phiSub: 1*n cell array of eigenfunctions evaluated at the t{i}'s
% sigma : variance of the measurement error
%         with measurement error assumption, it is > 0
%         without measurement error assumption, it is 0
% regular: 0 sparse and irregular data
%          2 completely balance data
% Output:
% invalid: 0 not valid
%          1 valid
% logLik1: marginal log-likelihood based on the inputs         
%
% See also convertMuPhi no_BIC no_AIC
function [invalid, logLik1] = getLogLik1(y,lambda, muSub, phiSub, sigma, regular)
   
     invalid = 0;
     LAMBDA = diag(lambda);
     logLik1 = 0;
     k = length(lambda);
     if regular == 2
        phiSub = phiSub(:,1:k);
        error0 = sigma*eye(length(muSub));
        Sigma_y = phiSub*LAMBDA*phiSub'+error0;
        detSigma_y = det(Sigma_y);
        if detSigma_y == 0
           invalid = 1;
           logLik1 = [];
           return;
        end
        for i = 1:length(y)           
           logLik1 = log(detSigma_y)+(y{i}-muSub)*pinv(Sigma_y)*(y{i}-muSub)' + logLik1;
        end      
     else
        for i = 1:length(y)
           phii = phiSub{i};
           phii = phii(:,1:k);
           mu_i = muSub{i};
           error0 = sigma*eye(length(mu_i));
           Sigma_y = phii*LAMBDA*phii'+error0;
           detSigma_y = det(Sigma_y);
           if detSigma_y == 0
             invalid = 1;
             logLik1 = [];
             return;
           else
             logLik1 = logLik1+log(detSigma_y)+(y{i}-mu_i)*pinv(Sigma_y)*(y{i}-mu_i)'; 
           end
        end
     end
end
