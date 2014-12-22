% ===========
% Description:
% ===========
% 
%            This function is used to compute pairwise distances for
%            sparsely observed functional trajectories. The output
%            distance is the estimated conditional expection of the L2
%            distance on observed data.
% 
% 
%            References: Yao, F., M\"{u}ller, H.G., Wang, J.L. (2005). Functional
%            data analysis for sparse longitudinal data.  J. American
%            Statistical Association 100, 577--590.
%            Peng, J., M\"{u}ller, H.G. (2008). Distance-based clustering of sparsely
%            observed stochastic processes, with applications to online auctions.
%            Annals of Applied Statistics 2, 1056--1077.
% 
%  
% ======
% Usage:
% ======
% 
% function [D1,D2,Y] = spadis(X)
% 
% ======
% Input: 
% ======
%      X:          a cell array containing the original data or the output
%                  of FPCA; if the original data are given, X{1} should be
%                  y, X{2} should be t and X{3} should be p as in the FPCA
%                  function; if user already used PACE, X should be a cell
%                  array containing all the outputs of the FPCA function.
% 
% =======
% Output:  
% =======  
%      D1:         n*n matrix where D1(i,j) is the estimated
%                  conditional expection of the L2 distance between subject
%                  i and j.  
%      D2:         n*n matrix where D2(i,j) is the Euclidean distance
%                  between estimated functional principal component socres
%                  of subject i and j.  
%      Y:          a cell array available when the original data are
%                  given, recording the output of PACE from the FPACE
%                  function.


function [D1,D2,Y] = spadis(X)

if length(X) == 3
    Y = FPCA(X{1},X{2},X{3});
else
    Y = X;
end
xi_var = getVal(Y,'xi_var');
xi_est = getVal(Y,'xi_est')';
if (size(xi_est,1) == 1)
  xi_est = [xi_est; zeros(1,size(xi_est,2))]; 
end
xi2 = sum(xi_est.*xi_est);
xi3 = xi_est'*xi_est; 
D2 = sqrt(repmat(xi2',[1,size(xi_est,2)]) + repmat(xi2,[size(xi_est,2),1]) - 2*xi3);
D2 = D2.*(1-eye(size(D2)));
tempVar = [];
for i = 1:length(xi_var)
    tempVar(i) = trace(xi_var{i});
end
D1 = sqrt(D2.^2+repmat(tempVar',[1,length(tempVar)])+repmat(tempVar,[length(tempVar),1]));
D1 = D1.*(1-eye(size(D1)));
if length(X)~=3 Y=[]; end

end
    
    
    
    
    
    
    
