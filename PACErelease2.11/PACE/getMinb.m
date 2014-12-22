%[dstar] = getMinb(t,out1,regular, npoly)
% This function approximates the minimum bandwidth choice for
% for the covariance function. 
% The starting bandwidth choice should be the minimum of this
% dstar and range(out1)/4.
% Input:
%     t: 1 * n cell array of time points for n subjects
%  out1: 1 * N vector of distinct unique time points from t
%  regular: 0 sparse and irragular
%           2 completely balance design
%  npoly: degree of polynomial to be fitted, by default, it is 1
% Output:
%  dstar: a positive value of the approximated minimum gap from the design plot 
%
% See also gcv_mullwlsn, cv_mullwlsn
function [dstar] = getMinb(t,out1,regular, npoly)

if nargin < 4
   npoly = 1;
end

%t1 = cputime;
if regular == 0
  [res] = designPlotCount(t,out1,1,1); %obtain the count matrix based on the design plot 
  %  fprintf(1,['Time after the count matrix: ' num2str(cputime-t1) '\n']);
  %  t1 = cputime;

  dstar = minb(out1,2+npoly);              %rough initial estimate of dstar based on 1-D sense
                                           %for at least 3 points
  res(:,[1 end]) = 1;
  res = res';
  ids = (res > 0);
  clear res;
  b = repmat(out1,length(out1),1)';
 
  dstar = max(dstar,max(diff(b(ids)))/2);

  clear b;
elseif regular == 1
  dstar = minb(out1,1+npoly)*2;
else
  %  t1 = cputime;
  dstar = minb(out1,2+npoly)*1.5;
end
%fprintf(1,['Time after the search of dstar: ' num2str(cputime-t1) '\n']);


end
