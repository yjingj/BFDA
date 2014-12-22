%Choose the best number of principal components through AIC method
function [no_opt,aic]=no_AIC(y,t,mu,bw_xcov,ngrid,regular,maxk,AB_method,method,shrink,out1,out21,kernel,error,cut,rcov,xcov,npoly)

if nargin < 18
  npoly = 1;
end
aic=Inf*ones(1,maxk);

%t1 = cputime;
if error == 1
   [invalid, sigma]=pc_covE(t,out1,bw_xcov,ngrid,cut,kernel, rcov, npoly);
else
   sigma = 0;
end
%t1 = cputime-t1;
%fprintf(1,['Time after pc_covE : ' num2str(t1) '\n']);

[lambda, phi, eigen, noeig]= getEigens(xcov,out1,out21,maxk);
[muSub, phiSub] = convertMuPhi(t,out1,mu,phi, regular);

if noeig < maxk
  fprintf(1, ['Warning: max number of PC for AIC selection is no more than ' num2str(noeig) ' ! You may want to increase maxk = ' num2str(maxk) ' to a larger number to include greater flexibility.\n']);
  maxk = noeig;
end

for k=1:maxk
%    t1 = cputime;
     sig = lambda(1:k);

%    fprintf(1,['Time after pc_est : ' num2str(t1) '\n']);
                                    
     if strcmp(AB_method,'full')
        [invalid, logLik] = getLogLik1(y,sig, muSub, phiSub, sigma,regular);               
        if invalid
	   fprintf(1, 'Warning: the covariance matrix of the estimated function is nearly singular! Reset to FVE method now!\n');
           no_opt = []; aic = [];
           return;
         end
     else
	 logLik=getLogLik2(y, t, sigma, k,method, shrink, regular, muSub, phiSub, sig);
     end
     aic(k)=logLik+2*k;  
     tmp = aic(~isinf(aic));
     if length(tmp) > 1 &&  k > 1
        if aic(k) > aic(k-1)
           break;
        end
     end       
end
aic = aic(~isinf(aic));
[temp,no_opt] = min(aic);

end
