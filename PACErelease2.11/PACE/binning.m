% This is a binning function.
% Observe data (x_i, y_i), i = 1, ..., n
% We want to convert the data into M bins or some bins with equidistant
% bin width to be h. The values of y_i that belonging to each bin with
% be replaced by their average.
%
% function [res] = binning(x, y, M, isMnumBin, nonEmptyOnly,a0,b0)
%
%======
%Input:
%======
% x        :  1 * n vector of values to be binned
% y        :  1 * n vector of values corresponding to x
% M        :  positive interger denotes the number of bins to be use
%             or
%             positive value denotes the width of each equidistant bin
%
% isMnumBin    :  1,  M denotes the number of bins
%                 0,  M denotes the width of each equidistant bin
%   
%  By default, isMnumBin is set to be 1, i.e., M is assumed to be number
%  of bins.  
%
% nonEmptyOnly:  1,  only non-empty bins will be output
%                0, all bins will be output
%
%  By default nonEmptyOnly is set to be 0, i.e., all bins
%  will be outputed.
%
%   a0    :   lower limit of domain of x 
%   b0    :   upper limit of domain of x
%   
%   By default, they are set to be min of x and max of x, respectively.
%======
%Output:
%======
% res : a struct object that contains midpoint, newy, count, numBin
%       and binWidth:
%        
%      i) midpoint: 1 * numBin vector of midpoints for each bin
%     ii) newy    : 1 * numBin vector of the corresponding y values
%                   for midpoint
%    iii) count   : 1 * numBin vector of counts in each bin
%     iv) numBin  : an integer denotes the number of bins
%      v) binWidth: a positive value denotes the width of each bin
%
% example:
% >>x = randn(100,1);             %generate 100 N(0,1) r.v.
% >>y = x+randn(100,1);           %generate another 100 N(0,1) r.v.
% >>M = 10;
% >>res = binning(x,y,M) % M denotes 10 bins includes all bins
% equivalently, 
% >>res = binning(x,y,M,1)
% 
% >> res = binning(x,y, 0.1,0,0)  % M denotes bin width
%                                 % include all bins
% equivalently,
% >> res = binning(x,y,0.1,0)   
% >> res = binning(x,y,0.1,0,1)   % M denotes bin width
%                                 % remove any empty bins
% >> res = binning(x,y,200,1,1)   % M denotes 200 bins
%                                 % remove any empty bins
%
% >> res = binning(x,y,10,1,1,0,5)% M is 10 bins
%                                 % remove any empty bins
%                                 % bin width is defined in
%                                 % [0,5] not [0,1]
function [res] = binning(x, y, M, isMnumBin, nonEmptyOnly,a0,b0)
                                                                
                                                                                             
  if M <= 0
     fprintf(1,'Error: M must be positive!\n');
     return;
  end
    
  if nargin == 2
     M = 10;
     isMnumBin = 1;
     a0 = min(x);
     b0 = max(x);                                                                                                                                          
  elseif nargin == 3
     M = ceil(M);
     isMnumBin = 1;
     nonEmptyOnly = 0;
     a0 = min(x);
     b0 = max(x);                                                                                                                                                
     if M == 1
         res = getResMisOne(x,y);
         return;
     end
  elseif nargin == 4
     nonEmptyOnly = 0;
     a0 = min(x);
     b0 = max(x);                                                                                                                                                           
     if isMnumBin 
       M = ceil(M);
       if M == 1
         res = getResMisOne(x,y);
         return;
       end
     end
  elseif nargin == 5
     a0 = min(x);
     b0 = max(x);                                                                                                                                                           
     if isMnumBin 
       M = ceil(M);
       if M == 1
         res = getResMisOne(x,y);
         return;
       end
    end
  elseif nargin == 6
      fprintf(1,'Error: b0 is missing!\n');
      return; 
  end

  if isMnumBin == 0
    h = M;
    if h >= range(x)
       res = getResMisOne(x,y,h);
       return;
    end
   
    xx = a0:h:b0;
    M = length(xx)-1;
  elseif isMnumBin

    if M == 1
         res = getResMisOne(x,y);
         return;
    end
    xx = linspace(a0,b0,M+1);
    h = xx(2)-xx(1);
  end

  N = length(xx);
  midpoint = xx(1:(N-1))+h/2;
  
  %xx(i-1) is like the lower limit
  %xx(i) is like the upper limit
  %in each bin, it includes the left end point 
  [newy, count] = getBins(x,y,xx,N);
  if nonEmptyOnly
    midpoint = midpoint(count > 0);
    newy = newy(count > 0);
    count = count(count>0);
    M = length(midpoint);
  end
  res = struct('midpoint', midpoint,'newy',newy, 'count', count,...
               'numBin', M, 'binWidth', h);

end


function [newy, count] = getBins(x,y, xx, N)

  count = zeros(1,N-1);
  newy = count;
 
   for i = 2:(N-1)
      ids = (x >= xx(i-1) & x < xx(i));
      if all(ids == 0)
       count(i-1) = 0;
       newy(i-1) = NaN;
     else
       count(i-1) = sum(ids);
       newy(i-1) = mean(y(ids));
     end
  end

  %for the last bin, it includes left and right end point
  ids = (x >= xx(i) & x <= xx(i+1));
  if all(ids == 0)
     count(i) = 0;
     newy(i) = NaN;
  else
     count(i) = sum(ids);
     newy(i) = mean(y(ids));
  end

end


function [res] = getResMisOne(x,y,h)
     r = range(x);
     if nargin == 2
       h = r;
     end   
     M = 1;
     midpoint = r/2;
     count = length(x);
     newy = mean(y);
     res = struct('midpoint', midpoint,'newy',newy, 'count', count,...
    'numBin', M, 'binWidth', h);
end

