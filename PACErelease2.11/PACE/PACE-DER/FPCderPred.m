%function [ypred, xi_new, xi_var] =  FPCderPred(yy, newy, newt, regular, idx)
%This function performs prediction for new y and t based on the returned
%fits from FPCder(). Unlike FPCderEval(), which is for the currently included subjects only,
%FPCderPred works for both new and currently included subjects.
% 
%Input yy:        an object that is returned by FPCder().
%
%Input newy:      1*m cell array of new measurements for new subjects
%
%Input newt:      1*m cell array of new time points for new subjects
%                 if all new subjects are evaluated at the same time,
%                 newt can be a row vector of time points 
%
%Input regular:   if the new time points are completely regular, then set regular = 2
%                 otherwise, set it to 0. 
%                 If it is [], use the same as defined in 'p'.  [Default] 
%                 When in doubt, set it to 0.    
%
%Input idx:       if not specified, predict all curves as specified in 'nder' [Default]
%                 or it can be a vector of 1,2,3 etc for predicting 
%                 curves correspond to nder(1),nder(2), nder(3) etc
%
%Output ypred:    1*length(idx) cell array of predicted measurements for new subjects
%                 ypred{j} are predicted curves for all new subjects with order nder(j)
%                 where j is the j-th element of the vector 'idx'.
%
%Output xi_new:   m x K matrix of new estimated FPC scores
%
%Output xi_var:   K*K matrix, Var(PC score)-Var(estimated PC score). The
%                 omega matrix in equation (7) of the paper, which is used
%                 to construct the point-wise C.I. for X_i(t)
%
% example:
% p = setDerOptions('nder',[0 1 2]);   %goal is to estimate curve and its 1st and 2nd derivative
% yy = FPCder(y,t,p);
% newy = {[1 2 3],4, [9 10]};
% newt = {[0.1 0.5 0.8], 1.1, [0.5 0.9]};
% [ypred, xi_new,xi_var] = FPCderPred(yy,p,newy,newt,[]); %use 'regular' defined in p
% ypred{1} : curve estimation for new subjects
% ypred{2} : 1st derivative estimation for new subjects
% ypred{3} : 2nd derivative estimation for new subjects
% 
% [ypred, xi_new,xi_var] = FPCderPred(yy,p,newy,newt,[],2); %use 'regular' defined in p
%                                                           %goal is estimating 2nd derivative
%  ypred{1} : 2nd derivative estimation for new subjects
% 
%  see also FPCder, FPCderEval
function [ypred, xi_new, xi_var] =  FPCderPred(yy, newy, newt, regular, idx)

  p = getVal(yy,'ops');
  if isempty(regular)
    regular = getVal(p,'regular');
  end  
  if nargin < 5 || idx < 1
     idx = 1:length(getVal(yy,'y_pred'));
  end

  if isnumeric(newt)
    if size(newt,1) > 1 && size(newt,2) == 1
       newt = newt';
    end
    newt = num2cell(repmat(newt,length(newy),1),2);
    regular = 2;
  end

  ypred = cell(1,length(idx));

  mu = getVal(yy,'mucopy');
  phi = getVal(yy,'phicopy');
  lambda = getVal(yy,'lambda');
  sigma = getVal(yy,'sigma');
  sigmanew = getVal(yy,'sigmanew');
  no_opt = getVal(yy,'no_opt');

  error = getVal(p,'error');
  method = getVal(p,'method');
  shrink = getVal(p,'shrink');
  out1 = getVal(yy,'out1copy');
  rho = getVal(yy,'rho_opt');
  n = length(getVal(yy,'y'));

  %obtain new scores
  [xi_new, xi_var,y0]=getXI(newy, newt, mu{1}, phi{1}, lambda, sigma, sigmanew, no_opt, error, method, shrink, out1, regular, rho);
  for j = 1:length(idx)
       if idx(j) == 1
           ypred{j} = y0;
       elseif idx(j) > 1
          [muSub, phiSub] = convertMuPhi(newt, out1, mu{idx(j)}, phi{idx(j)}, regular); 
          if regular == 2
              tmp = repmat(muSub, n,1)+xi_new*phiSub';
              ypred{j} = num2cell(tmp,2);
           else
              ypredj = cell(1,n);
              for i = 1:n
		ypredj{i} = muSub{i}+ xi_new(i,:)*phiSub{i}';          
              end
	      ypred{j} = ypredj;
           end
       end
  end
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
