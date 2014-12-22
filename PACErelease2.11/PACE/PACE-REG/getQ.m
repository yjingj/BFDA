% function [Q] = getQ(fitted_y,yy , y, t_y,isYFun)
% Input:
%  fitted_y: a cell array of the same length as y. fitted_y{i} is the
%            fitted values of the i-th subject (at the original observation time points).
%        yy: an object that contains the returned values from FPCA(y,t_y,param_Y) or "[]" is y is scalar response
%            see PCA() for more details, use getVal() to retrieve the individual values.
%         y: the original input response
%       t_y: the original input observation time for response  
%    isYFun: a logical value, isYFun = 1 (Y is a functional response)
%                             isYFun = 0 (Y is a scalar response)
%     dense: a logical value for functional Y
%            dense = 1, Y is dense; dense = 0, Y is sparse.
% Output:
%    Q:   Quasi R-square; a measure of variation explained in the response.
%         If Y is a scalar, 
%                Q = 1-sum((Y_i-Yhat_i)^2)/sum((Y_i-mean(Y))^2)
%         If Y is a function,
%                Q = 1-sum((Y_i-Yhat_i)^2)/n_i)/sum((Y_i-mean(Y))^2/n_i),   if dense = 0
%                Q = 1-\int((Y_i-Yhat_i)^2)/\int((Y_i-mean(Y))^2),          if dense = 1

function [Q] = getQ(fitted_y, yy, y, t_y, isYFun, dense)

if nargin < 6
    dense = 0;
end    
    
if isYFun == 1

    muy = getVal(yy,'mu');
    out1y = getVal(yy,'out1');
    ncohort = length(y);
    v1 = zeros(1,ncohort);
    v2 = zeros(1,ncohort);    
    if dense == 0

        for i = 1:ncohort
            v1(i) = sum((y{i}-fitted_y{i}).^2)/length(y{i});
            v2(i) = sum((y{i}-interp1(out1y,muy,t_y{i},'spline')).^2)/length(y{i});
        end
        
    else
        
        for i = 1:ncohort
            v1(i) = trapz(t_y{i},(y{i}-fitted_y{i}).^2);
            v2(i) = trapz(t_y{i},(y{i}-interp1(out1y,muy,t_y{i},'spline')).^2);
        end
        
    end
    
    Q = 1-sum(v1)/sum(v2);

else

    Q = 1-sum((y-fitted_y).^2)/sum((y-mean(y)).^2);

end

if Q < 0
    fprintf(1,'Warning: Negative Quasi R-square!\n');
end


