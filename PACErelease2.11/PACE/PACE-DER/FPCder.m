%[X] = FPCder(y,t,p)
%
%===========
%Description:      
%=========== 
%
%            This is the main program to perform Functional Principal Component
%            Analysis (FPCA) via PACE for estimation of curves and their derivatives. 
%            The principal component scores can be estimated through conditional 
%            expectation or via classical integration. For the latter, one can 
%            choose a shrinkage method for estimated scores.
%
%======
%Input: 
%======
%      y:          1*n cell array, y{i} is the vector of measurements for the ith subject,
%                  i=1,...,n.
%      t:          1*n cell array, t{i} is the vector of time points for the ith subject on which
%                  corresponding measurements y{i} are taken, i=1,...,n.
%      p:          a struct obtained from setOptions.m sets the rest of arguments for PCA.m
%                  ex:
%                  >> p = setDerOptions();
%=======
%Output:  
%=======  
%      X:         a cell array that contains all returned values from PCA.m
%                 where the last element contains the names of X
% 
%   To see the names for X, type names(X)
%   To get individual value back, type getVal(varname, X)
%   To see an example, check with exampleDer.m
%   See also PCder, setDerOptions, showDerOptionNames, exampleDer, names
%
function [X] = FPCder(y,t,p)

  if nargin == 2
    p = setDerOptions();
  end
 %addpath('./PACE');
  if isempty(p.yname)
     p.yname = inputname(1);
  end
[no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov,xcovfit,AIC,BIC,...
 FVE,y_pred,y_predOrig,y_predDense,out1,out21,y,t, rcov, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops]...
=PCder(y,t,p.bwmu, p.bwmu_gcv, p.bwxcov, p.bwxcov_gcv, ...
p.ntest1,p.ngrid1,p.selection_k, p.FVE_threshold, p.maxk,...
p.control,p.regular,p.error, p.ngrid,p.method,p.shrink,p.newdata, p.kernel, ...
p.npoly, p.nder, p.method_int, p.numBins, p.yname, p.screePlot, p.designPlot, p.corrPlot, p.rho, p.verbose);



Xnames = {'no_opt','sigma','lambda','phi','eigen','xi_est','xi_var','mu','muDense','bw_mu','xcov','bw_xcov',...
	  'xcovfit','AIC','BIC','FVE','y_pred','y_predOrig','y_predDense','out1','out21', 'y','t','rcov', 'regular', 'rho_opt','sigmanew','mucopy',...
	  'phicopy','eigencopy','out1copy', 'out21copy', 'xcovcopy', 'xcovfitcopy','xcorr','ops','names'};
X = {no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov,...
     xcovfit,AIC,BIC,FVE,y_pred,y_predOrig,y_predDense,out1,out21, y,t,rcov, regular, rho_opt, sigmanew, mucopy,...
     phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops, Xnames};  

end


