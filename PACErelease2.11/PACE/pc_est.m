function [xi_est, xi_var, lambda, phi, eigen, noeig, xcovfit,y_predOrig, rho_opt, sig1]=pc_est(y, t, mu, xcov, sigma, noeig, error, method, shrink, out1, out21,regular, rho, verbose)

[lambda, phi, eigen, noeig] = getEigens(xcov,out1,out21,noeig,1);
%calculate the fitted covariance 
%based on lambda and eigen
%fprintf(1,['Total final selected number of PC: ' num2str(noeig) ' .\n']);
xcovfit = zeros(size(eigen,1),size(eigen,1));
for i = 1:noeig
   xcovfit = xcovfit + lambda(i)* eigen(:,i)*eigen(:,i)';   
end

[xi_est, xi_var,y_predOrig, rho_opt, sig1]=getScores(y, t, mu, phi, lambda, sigma, noeig, error, method, shrink, out1, regular, rho, verbose);
