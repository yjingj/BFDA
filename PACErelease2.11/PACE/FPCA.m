%[X] = FPCA(y,t,p)
%This function calls PCA.m.
%======
%Input: 
%======
%      y:          1*n cell array, y{i} is the vector of measurements for the ith subject,
%                  i=1,...,n.
%      t:          1*n cell array, t{i} is the vector of time points for the ith subject on which
%                  corresponding measurements y{i} are taken, i=1,...,n.
%      p:          a struct obtained from setOptions.m sets the rest of arguments for PCA.m
%                  ex:
%                  >> p = setOptions();
%=======
%Output:  
%=======  
%      X:         a cell array that contains all returned values from PCA.m
%                 where the last element contains the names of X
% 
%   To see the names for X, type names(X)
%   To get individual value back, type getVal(X,varname)
%   To see an example, check with example.m
%   See also PCA, example, names
%
function [X] = FPCA(y,t,p)

  if nargin == 2
    p = setOptions();
  end
 %addpath('./PACE');
  if isempty(p.yname)
     p.yname = inputname(1);
  end
  [no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov,xcovfit,AIC,BIC,FVE,y_pred,y_predOrig,y_predDense,out1,...
	 out21,y,t, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops]...
=PCA(y,t,p.bwmu, p.bwmu_gcv,p.bwxcov, p.bwxcov_gcv,p.ntest1,p.ngrid1,p.selection_k, p.FVE_threshold, p.maxk,...
     p.control,p.regular,p.error,p.ngrid,p.method,p.shrink,p.newdata, p.kernel, p.numBins, p.yname, p.screePlot, p.designPlot, p.corrPlot, p.rho, p.verbose, p.xcov);

Xnames = {'no_opt','sigma','lambda','phi','eigen','xi_est','xi_var','mu','muDense','bw_mu','xcov','bw_xcov','xcovfit',...
	  'AIC','BIC','FVE','y_pred','y_predOrig','y_predDense','out1','out21', 'y','t','regular','rho_opt', 'sigmanew'...
	  'mucopy', 'phicopy','eigencopy','out1copy', 'out21copy', 'xcovcopy', 'xcovfitcopy','xcorr','ops', 'names'};
X = {no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov,xcovfit,AIC,BIC,FVE,...
     y_pred,y_predOrig,y_predDense,out1,out21, y,t, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops,Xnames};  
end


