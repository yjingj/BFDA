%[y] = mysample(x,n,withReplace)
%sample n objects from x with or without replacement
%x:           a vector or cell array to be sampled
%n:           sample n elements from x
%             when sampling with replacement, n is allow to be larger than
%             the size of x;
%             when sampling without replacement, n is at most the size of x
%withReplace: 1 sample x with replacement
%             0 sample x without replacement
%
%Example: 
%1) x = randn(10,1);
%   y = mysample(x,20,1);
%2) y = mysample(x,6,0);
%3) x = {1,1:10, 'something',[1 3; 5 4]}
%   y = mysample(x);
function [y] = mysample(x,n, withReplace)
     if nargin == 1
        n = length(x);
        withReplace = 1;
     elseif nargin < 3
        withReplace = 1;
     end
     rand('twister',sum(10000*clock));
     if n <= 0 
        error('n must be a positive integer!');
     elseif n > length(x) && withReplace == 0
        error('When sample without replacement, n must be no larger than x!');
     end
     [something,ind] = sort(rand(length(x),1));
     if withReplace == 0
         y = x(ind(1:n));
     else
         ids = ceil(length(x)*rand(n,1));
         y = x(ids(1:n));
     end

end
