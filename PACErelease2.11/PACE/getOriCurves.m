function [y_predOrig,xi_est]=getOriCurves(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular)

  [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular);
  ncohort = length(y);
  LAMBDA = diag(lambda);

  if error==1 
      
      sigma1 = sig1;

      if regular == 2 && strcmp(method,'CE')
 
         yy= reshape(cell2mat(y), length(y{1}), ncohort)';
 
         error0 = sigma1*eye(length(t{1}));
         A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub'+error0);
         MU = repmat(muSub, ncohort,1);
         B = yy-MU;
         xi_est = (A*B')';
         y_predOrig = MU+xi_est*phiSub';
         y_predOrig = num2cell(y_predOrig,2);

      else
	     y_predOrig = cell(1,ncohort);
	     xi_est = zeros(ncohort, noeig);
	     zeta_est = xi_est;
	     phii= phiSub;
	     mu_i = muSub;
	     for i = 1:ncohort

            if regular ~= 2
	           phii = phiSub{i};
	           mu_i = muSub{i};
            end
	        yi= y{i};
	        if strcmp(method,'CE')              
                  error0=sigma1*eye(length(yi));
		  A = LAMBDA*phii'*pinv(phii*LAMBDA*phii'+error0);
                  xi_est(i,:)=(A*(yi-mu_i)')';

	        elseif strcmp(method,'IN')
                  m=length(yi); 
                  for k=1:noeig         
                    prod=(yi-mu_i).*phii(:,k)';
                    if shrink == 0
                       %xi_est(i,k) = romb(t{i},prod);
                       xi_est(i,k) = trapz(t{i},prod);
                    else
                       %zeta_est(i,k) = romb(t{i},prod);
                       zeta_est(i,k) = trapz(t{i},prod);
                       xi_est(i,k)=lambda(k)*zeta_est(i,k)/(lambda(k)+sigma/m);
                    end
                  end
                end
                y_predOrig{i} = mu_i+xi_est(i,:)*phii';
         end
      end
  elseif error==0 
      if regular == 2 && strcmp(method ,'CE')
          yy= reshape(cell2mat(y), length(y{1}), ncohort)';
          A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub');
          MU = repmat(muSub, ncohort,1);
          B = yy-MU;
          xi_est = (A*B')';
          y_predOrig = MU+xi_est*phiSub';
          y_predOrig = num2cell(y_predOrig,2);
      else
          y_predOrig = cell(1,ncohort);
       	  xi_est = zeros(ncohort, noeig);
	      phii= phiSub;
	      mu_i = muSub;
	      for i = 1:ncohort
             if regular ~= 2
		        phii = phiSub{i};
		        mu_i = muSub{i};
             end
	         yi= y{i};
	         if strcmp(method,'CE')
			 A = LAMBDA*phii'*pinv(phii*LAMBDA*phii');
		         xi_est(i,:)=(A*(yi-mu_i)')';
	         elseif strcmp(method,'IN')
		         for k=1:noeig
		            prod=(yi-mu_i).*phii(:,k)';
                            %xi_est(i,k) = romb(t{i},prod);
                            xi_est(i,k) = trapz(t{i},prod);
                         end
                 end
                 y_predOrig{i} = mu_i+xi_est(i,:)*phii';
          end
      end
  end
end
