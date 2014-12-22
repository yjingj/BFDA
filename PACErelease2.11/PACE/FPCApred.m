%function [ypred, xi_new, xi_var] =  FPCApred(yy, newy, newt, regular)
%This function performs prediction for new y and newt based on the returned
%fits from FPCA(). Unlike FPCAeval, which is for the currently included subjects only,
%FPCApred works for both new and currently included subjects.
% 
%Input yy:        an object that is returned by FPCA().
%Input newy:      1*m cell array of new measurements for new subjects
%Input newt:      1*m cell array of new time points for new subjects
%                 if all new subjects are evaluated at the same time,
%                 newt can be a row vector of time points for one subject
%Input regular:   0,1,2 (see PCA() for more details), if in doubts, set regular = 0
%                 if not specified, then, the program uses the same 'regular' as specified in 'p'.
%Output ypred:    1*m predicted measurements for new subjects
%Output xi_new:   m x K matrix of new estimated FPC scores
%Output xi_var:   K*K matrix, Var(PC score)-Var(estimated PC score). The
%                 omega matrix in equation (7) of the paper, which is used
%                 to construct the point-wise C.I. for X_i(t)
% example:
% p = setOptions();
% yy = FPCA(y,t,p);
% newy = {[1 2 3],4, [9 10]};
% newt = {[0.1 0.5 0.8], 1.1, [0.5 0.9]};
% [ypred, xi_new,xi_var] = FPCApred(yy,newy,newt); %use 'regular' defined in p
%
%  see also FPCA, FPCAeval
function [ypred, xi_new, xi_var] =  FPCApred(yy,newy, newt, regular)

  p = getVal(yy,'ops');
  if nargin < 4
    regular = getVal(p,'regular');
  end
  
  if isnumeric(newt)
    if size(newt,1) > 1 && size(newt,2) == 1
       newt = newt';
    end
    newt = num2cell(repmat(newt,length(newy),1),2);
    regular = 2;
  end

  mu = getVal(yy,'mucopy');
  phi = getVal(yy,'phicopy');
  phi = phi(:,1:getVal(yy,'no_opt'));
  lambda = getVal(yy,'lambda');
  lambda = lambda(1:getVal(yy,'no_opt'));
  sigma = getVal(yy,'sigma');
  sigmanew = getVal(yy,'sigmanew');
  no_opt = getVal(yy,'no_opt');

  error = getVal(p,'error');
  method = getVal(p,'method');
  shrink = getVal(p,'shrink');
  out1 = getVal(yy,'out1copy');
  rho = getVal(yy,'rho_opt');

  [xi_new, xi_var,ypred]=getXI(newy, newt, mu, phi, lambda, sigma, sigmanew, no_opt, error, method, shrink, out1, regular, rho);

end

function [xi_est, xi_var,y_predOrig]=getXI(y, t, mu, phi, lambda, sigma, sigma_new, noeig, error, method, shrink, out1, regular, rho)

  [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular);
  ncohort = length(y);
  LAMBDA = diag(lambda);
  if strcmp(method,'IN')
     xi_var = [];
  else
     xi_var = cell(1,ncohort);
  end
  %fprintf(1,'Start calculating the PC scores \n');
  %update \xi 

  if error==1 
      
      if sigma_new < rho
	   sigma1 = rho;
      else
	   sigma1 = sigma_new;
      end

      if regular == 2 && strcmp(method,'CE')
 
         yy= reshape(cell2mat(y), length(y{1}), ncohort)';
 
         error0 = sigma1*eye(length(t{1}));
         A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub'+error0);
         MU = repmat(muSub, ncohort,1);
         B = yy-MU;
         xi_est = (A*B')';
         y_predOrig = MU+xi_est*phiSub';
         y_predOrig = num2cell(y_predOrig,2);
         C = LAMBDA-A*(LAMBDA*phiSub')';
         for i = 1:ncohort
            xi_var{i} = C;
         end     

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
                  xi_var{i}=LAMBDA-A*(LAMBDA*phii')';
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
          C = LAMBDA-A*(LAMBDA*phiSub')';
          for i = 1:ncohort
    	     xi_var{i} = C;
          end
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
		         xi_var{i}=LAMBDA-A*(LAMBDA*phii')';
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


function [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular)

     if regular == 2
	muSub = interp1(out1,mu,t{1},'spline');
	phiSub = interp1(out1,phi,t{1},'spline');
        if size(phi,2) == 1
           phiSub = phiSub';
        end
     else
        muSub = cell(1,length(t));
        phiSub = cell(1,length(t));
        for i = 1:length(t)
	    muSub{i} = interp1(out1,mu,t{i},'spline');
	    phiSub{i} = interp1(out1,phi,t{i},'spline');
            if size(phi,2) == 1
               phiSub{i} = phiSub{i}';
            end
        end
     end
end
