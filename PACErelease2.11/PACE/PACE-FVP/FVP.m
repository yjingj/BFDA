% ============
% Description:
% ============
% 
%   This is the main function to perform Functional Variance Process Analysis
% 
% Includes the following steps:
% 
% 1) PACE Estimation for the original processes Y.
% 
% 2) Calculation of the raw residual process from step 1.
% 
% 3) PACE Estimation for the functional variance processes.
%
% See reference:  M\"uller, H.G., Stadtm\"uller, U., Yao, F. (2006). Functional variance processes. 
%                 Journal of the American Statistical Association 101, 1007-1018. 
% ========
% Usage:
% ========
% 
%     [rr,sigma,yy] = FVP(y,t,p)
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
%                  >> p = setOptions();
%=======
%Output:  
%=======  
%     rr:          a cell array that contains all returned values from PCA.m for the functional 
%                  variance processes.
%  sigma:          the variance estimator of the functional variance process.
%     yy:          a cell array that contains all returned values from PCA.m for the original 
%                  functional processes.
%      
%
%   To see the names for rr, type names(rr)
%   To get individual value back, type getVal(rr,varname)
%   To see an example, check with example.m
%   See also FPCA, PCA, example, names
%
% Note: Squared residuals are adjusted by adding the 10th percentile of % the pooled sample before taking log.
%
function [rr,sigma,yy] = FVP(y,t,p)
    q=.1;
    yy = FPCA(y,t,p);
    y_predOrig = getVal(yy,'y_predOrig');
    r = cell(1,length(y));
    n = length(y_predOrig);
    for i=1:n
        r{i} = (y{i}-y_predOrig{i}).^2;
    end
    delta = quantile(cell2mat(r),q);
    for i=1:n
        r{i} = log((y{i}-y_predOrig{i}).^2+delta);
    end
    rr = FPCA(r,t,p);
    sigma = getVal(rr,'sigma');
end