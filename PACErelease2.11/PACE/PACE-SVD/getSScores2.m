%This function is used by cv_srho.m and cv_sbw.m only
%It calculates the scores for leaving one random t_ij out for a given rho.
%It only applies to the case where error = 1
function cv = getSScores2(res, x, t_x,  y, t_y, nsvd, muxSub, scxSub, muySub, scySub, sig1, method, shrink, regular, rho, subID, tjID)
 
%   mu_x=getVal(getVal(res,'xx'),'mucopy');
%   sc_x=getVal(res,'sc_x');
  lambda_x = getVal(res,'lambda_x');
%   out1x = getVal(res,'out_x');
  sigma_x = getVal(getVal(res,'xx'),'sigma');
%   mu_y=getVal(getVal(res,'yy'),'mucopy');
%   sc_y=getVal(res,'sc_y');
  lambda_y = getVal(res,'lambda_y');
%   out1y = getVal(res,'out_y');
  sigma_y = getVal(getVal(res,'yy'),'sigma');
  lambda = getVal(res,'lambda');
  
  cv = 0;
  ncohort = length(subID); %here ncohort is the number of subjects with ni >= 2 only
  tjID = tjID(subID,:);
%  if error==1 
      
      sigma1 = max(rho, sig1);

      %if regular == 2 && strcmp(method,'CE')
 
      %   yy= reshape(cell2mat(y), length(y{1}), ncohort)';
 
      %   error0 = sigma1*eye(length(t{1}));
      %   A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub'+error0);
      %   MU = repmat(muSub, ncohort,1);
      %   B = yy-MU;
      %   xi_est = (A*B')';
      %   y_predOrig = MU+xi_est*phiSub';
      %   y_predOrig = num2cell(y_predOrig,2);
 
      %else
	     %y_predOrig = cell(1,ncohort);
	     xi_est = zeros(ncohort, nsvd);
         yi_est = xi_est;
         zeta_xest = xi_est;
         zeta_yest = xi_est;
% 	     scxi= scxSub;
% 	     mux_i = muxSub;
% 	     scyi= scySub;
% 	     muy_i = muySub;

         for i = 1:ncohort

             if regular ~= 2
                 scxi = scxSub{subID(i)};
                 mux_i = muxSub{subID(i)};
                 scyi = scySub{subID(i)};
                 muy_i = muySub{subID(i)};
             else
                 scxi = scxSub;
                 mux_i = muxSub;
                 scyi = scySub;
                 muy_i = muySub;
             end

             xi= x{subID(i)};
             xij = xi(tjID(i,1));
             xi(tjID(i,1)) = [];
             txi = t_x{subID(i)};
%              txij = txi(tjID(i,1));
             txi(tjID(i,1)) = [];
             mux_ij = mux_i(tjID(i,1));
             scxi_j = scxi(tjID(i,1),:);
             mux_i(tjID(i,1)) = [];
             scxi(tjID(i,1),:) = [];
             yi= y{subID(i)};
             yij = yi(tjID(i,2));
             yi(tjID(i,2)) = [];
             tyi = t_y{subID(i)};
%              tyij = tyi(tjID(i,2));
             tyi(tjID(i,2)) = [];
             muy_ij = muy_i(tjID(i,2));
             scyi_j = scyi(tjID(i,2),:);
             muy_i(tjID(i,2)) = [];
             scyi(tjID(i,2),:) = [];

             if strcmp(method,'CE')

                 A22 = scxi*diag(lambda_x)*scxi'+sigma1(1)*eye(length(xi));
                 B22 = scxi*diag(lambda)*scyi';
                 D22 = scyi*diag(lambda_y)*scyi'+sigma1(2)*eye(length(yi));
                 A22i = pinv(A22-B22*pinv(D22)*B22');
                 D22i = pinv(D22-B22'*pinv(A22)*B22);
                 B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
                 sig22i = [A22i B22i; B22i' D22i];
                 sig12 = [diag(lambda_x)*scxi' diag(lambda)*scyi'; diag(lambda)*scxi' diag(lambda_y)*scyi'];
                 mu_ss = sig12 * sig22i * ([xi yi] - [mux_i muy_i])';
                 xi_est(i,:) = mu_ss(1:nsvd);
                 yi_est(i,:) = mu_ss((nsvd+1):end);

             elseif strcmp(method,'IN')

                 mx=length(xi);
                 my=length(xi);

                 for k=1:nsvd
                     prodx=(xi-mux_i).*scxi(:,k)';
                     prody=(yi-muy_i).*scyi(:,k)';
                     if shrink == 0
                         %xi_est(i,k) = romb(ti,prod);
                         xi_est(i,k) = trapz(txi,prodx);
                         yi_est(i,k) = trapz(tyi,prody);
                         %zeta_est(i,k)=trapzoid(prod,ti);
                     else
                         %zeta_est(i,k) = romb(t{i},prod);
                         zeta_xest(i,k) = trapz(txi,prodx);
                         xi_est(i,k)=lambda_x(k)*zeta_xest(i,k)/(lambda_x(k)+sigma_x/mx);
                         zeta_yest(i,k) = trapz(tyi,prody);
                         yi_est(i,k)=lambda_y(k)*zeta_yest(i,k)/(lambda_y(k)+sigma_y/my);
                     end                     
                 end
             end
             %y_predOrig{i} = mu_i+xi_est(i,:)*phii';
             cv = cv+(mux_ij+xi_est(i,:)*scxi_j' - xij)^2+(muy_ij+yi_est(i,:)*scyi_j' - yij)^2;
         end
      %end
% elseif error==0 
%     if regular == 2 && strcmp(method ,'CE')
%        yy= reshape(cell2mat(y), length(y{1}), ncohort)';
%         A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub');
%         MU = repmat(muSub, ncohort,1);
%         B = yy-MU;
%         xi_est = (A*B')';
%         y_predOrig = MU+xi_est*phiSub';
%         y_predOrig = num2cell(y_predOrig,2);
%
%     else
%         y_predOrig = cell(1,ncohort);
%      	  xi_est = zeros(ncohort, noeig);
%	      phii= phiSub;
%	      mu_i = muSub;
%	      for i = 1:ncohort
%                  if regular ~= 2
%		     phii = phiSub{i};
%		     mu_i = muSub{i};
%                  end
%	          yi= y{i};
%	          if strcmp(method,'CE')
%			 A = LAMBDA*phii'*pinv(phii*LAMBDA*phii');
%		         xi_est(i,:)=(A*(yi-mu_i)')';
%      	          elseif strcmp(method,'IN')
%		         for k=1:noeig
%		            prod=(yi-mu_i).*phii(:,k)';
%                            %xi_est(i,k)=trapzoid(prod,t{i});
%                            %xi_est(i,k) = romb(t{i},prod);
%                            xi_est(i,k) = trapz(t{i},prod);
%                         end
%                  end
%                  y_predOrig{i} = mu_i+xi_est(i,:)*phii';
%              end
%      end
% end

end
