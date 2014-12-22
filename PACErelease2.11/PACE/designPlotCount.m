%[res] = designPlotCount(t, out1, noDiagonal, isIndicator)
% This function is used in gcv_mullwlsn.m to create count
% matrix based on observed pairs of time points for the raw
% covariance.
% Input t:  n * 1 array contains time points for n subjects
% Input out1: 1 * N vector contains sorted unique time points from t
% Input noDiagonal: 1 set diagonal count as 0
%                   0 don't set diagonal count as 0
% Input isIndicator: 1 the resulting matrix has 1 indicates there exists points for out1(i) and out1(j)
%                    0 the resulting matrix contains counts of points for out1(i) and out1(j)
% Output res: N * N matrix contains count for each distinct pairs of 
%             time points
function [res] = designPlotCount(t,out1, noDiagonal,isIndicator)
N = length(out1); %number of unique time points
%res = zeros(N,N);
res = sparse(N,N);
for i = 1:length(t)
     cur = t{i};   %time points from subject i
     curidx = searchID(cur,out1);
     if isIndicator == 1
          res(curidx,curidx) = 1;
     else
          res(curidx,curidx) = res(curidx,curidx)+1;
     end
end
if noDiagonal == 1
   tmp = diag(res);
   res = res-diag(tmp);   
end
end 
function [id] = searchID(cur,out1)
   
   n = length(cur);
   id = zeros(n,1);
   for i = 1:n
      id(i) = find(out1 == cur(i));
   end
end

