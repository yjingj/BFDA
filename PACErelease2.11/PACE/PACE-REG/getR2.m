% function [r2] = getR2(b, xx, yy, y, K_x, K_y, isYFun)
% Input:
%    b:    a K_x*K_y matrix; estimated coefficients of beta in the basis represatation 
%          returned value from getBeta(...)
%    xx:   an object that contains the returned values from FPCA(x,t_x,param_X)
%          see PCA() for more details, use getVal() to retrieve the individual values.
%    yy:   an object that contains the returned values from FPCA(y,t_y,param_Y) or "[]" is y is scalar response
%          see PCA() for more details, use getVal() to retrieve the individual values.
%     y:   the original input response; used when calculating prediction
%          R-square for scalar case
%   K_x:   positive integer; number of principal components of predictor x used in regression
%          must be smaller than or equal to the maximum number of principal components give by FPAC
%   K_y:   positive integer; number of principal components of response y used in regression
%          must be smaller than or equal to the maximum number of principal components give by FPAC
% isYFun:  a logical value, isYFun = 1 (Y is a functional response)
%                           isYFun = 0 (Y is a scalar response)
% 
% Output:
%    r2:   functional R-square

function [r2] = getR2(b, xx, yy, y, K_x, K_y, isYFun)

    lambda_x = getVal(xx,'lambda'); 
    lambda_X = lambda_x(1:K_x);

    if isYFun ==1
        
        lambda_y = getVal(yy,'lambda'); 
        lambda_Y = lambda_y(1:K_y);    
        r2 = sum(lambda_X*(b.^2))/sum(lambda_Y);
        
    else
 
        r2 = sum((b.^2).*lambda_X)/var(y);
        
    end
    
    if r2>0.99
        fprintf(1,'Warning: Too large functional R^2. The model may be overfitted!\n');
        r2=0.99;
    end

end
