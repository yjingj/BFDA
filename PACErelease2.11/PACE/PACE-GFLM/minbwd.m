% function used to find the minimum bandwidth choice where local window
% contains at least "numPoints" points 
% Input x  : n x 1 vector
% Input numPoints: an integer specifying the number of points in a local window
% for local weighted constant, numPoints is at least 1
% for local weighted linear, numPoints is at least 2
% Output b: the minimum bandwidth choice for vector x

function [b] = minbwd(x,numPoints)

n = length(x);
x = sort(x);
if numPoints > 1
    b = max(x(numPoints:n)-x(1:(n-numPoints+1)));
else
    b = max((x(2:n)-x(1:(n-1)))/2);
end

end