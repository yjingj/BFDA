%This function obtain the conditional log-likelihood for AIC1 or BIC1
%It is used in no_AIC or no_BIC functions.
%Only use with measurement error assumption
%
% [logLik2]=getLogLik2(y, t, sigma, noeig,method, shrink, regular, muSub, phiSub, lambda)
% Input:
%     y: 1*n cell array of observed repeated measurements for n subjects
%     t: 1*n cell array of observed time points for n subjects
% sigma : variance of the measurement error, positive value
% noeig: positive integer of No. of PC
% method: 'CE' or 'IN' method to estimate the scores
%          conditional expectation or integration method
% shrink: 1 shrinkage estimate
%         0 no shrinkage
% regular: 0 sparse and irregular data
%            2 completely balance data
%  muSub: 1*n cell array of mean functions evaluated at the t{i}'s
% phiSub: 1*n cell array of eigenfunctions evaluated at the t{i}'s
% lambda: 1*k vector of eigenvalues
% Output:
% logLik2: conditional log-likelihood based on the inputs
%
% See also convertMuPhi no_BIC no_AIC

function [logLik2]=getLogLik2(y, t, sigma, noeig,method, shrink, regular, muSub, phiSub, lambda)

%fprintf(1,'Start calculating the PC scores \n');
%update \xi 
LAMBDA = diag(lambda);
k = length(lambda);

if regular == 2 && strcmp(method,'CE')
      phiSub = phiSub(:,1:k);
      yy= reshape(cell2mat(y), length(y{1}), length(y))';
      error0 = sigma*eye(length(t{1}));
      A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub'+error0);
      MU = repmat(muSub, length(y),1);
      B = yy-MU;
      xi_est = (A*B')';
      y_pred = MU+xi_est*phiSub';
      res = yy-y_pred;
      res = (res.^2)/sigma;
      logLik2 = sum(res(:));

else

      xi_est=zeros(1,noeig);
      zeta_est=xi_est;
      logLik2 = 0;
      phii = phiSub;
      mu_i = muSub;
      for i=1:length(y)
         if regular ~= 2
            phii = phiSub{i};
            mu_i = muSub{i};
         end
         phii = phii(:,1:k);
         yi = y{i};
         if strcmp(method,'CE')
            error0=sigma*eye(length(yi));
            A = LAMBDA*phii'*pinv(phii*LAMBDA*phii'+error0);
            xi_est=(A*(yi-mu_i)')'; 
         elseif strcmp(method, 'IN')
            m=length(yi); 
            for k=1:noeig         
                prod=(yi-mu_i).*phii(:,k)';
                zeta_est(k)=trapz(t{i},prod);
                xi_est(k)=lambda(k)*zeta_est(k)/(lambda(k)+sigma/m);
                if shrink==0
                    xi_est(k)=zeta_est(k);
                end
            end
         end
         y_pred = mu_i+xi_est*phii';
         logLik2 = (yi-y_pred)*(yi-y_pred)'/sigma+logLik2;        
      end

end

end

