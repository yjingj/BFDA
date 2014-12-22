%[isYFun] = checkY(y)
% This function checks whether response variable y is a 
% function (cell array) or a scalar response (vector).
% Input: y is a cell array or a vector
% Output: isYFun  a logical value, isYFun = 1 (Y is a functional response)
%                                  isYFun = 0 (Y is a scalar response)
function [isYFun] = checkY(y)

     if iscell(y) == 1
        isYFun = 1;
     elseif isnumeric(y) == 1
        isYFun = 0;
     else
        fprintf(1,'Error: y must be either a cell array or a numeric vector!\n');
     end

end
