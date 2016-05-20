function [sigma] = topdm(sig) 

% FUNCTION THAT TRANSFORMS A NON POSITIVE DEFINITE SYMMETRIC
% MATRIX -E.G. A COVARIANCE MATRIX- TO A 
% POSITIVE DEFINITE SYMMETRIC MATRIX. 


EPS =  (10^-6);      %SET THE VALUE TO PLACE INSTEAD OF ZERO OR NEGATIVE  
                  %EIGENVALUES  
ZERO = (10^-7);   %SET THE VALUE TO LOOK FOR 

sigma = sig; 
[~, err] = cholcov(sigma, 0); %ERR<>0 MEANS SOME EIGENVALUES <=0 
                              %So... SIGMA IS NOT POSITIVE DEFINITE 
                              %So... SIGMA IS NOT INVERTIBLE 

if (err ~= 0) 
    [v d] = eig(sigma); %CALCULATE EIGENVECTOR AND EIGENVALUES 
    d=diag(d);          %GET A VECTOR OF EIGENVALUES 
    d(d<ZERO)=EPS; %FIND ALL EIGENVALUES<=ZERO AND CHANGE THEM FOR EPS 
    d=diag(d);      %CONVERT VECTOR d INTO A MATRIX 
    sigma = v*d*v'; %RECOMPOSE SIGMA MATRIX USING EIGENDECOMPOSITION 
                    %WHY? SIGMA IS SIMETRIC AND V IS ORTHONORMAL 
end 

