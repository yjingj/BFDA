%Select best number of principal components based on fraction of
%variance explained (FVE).
%[no_opt, FVE, d] = no_FVE(xcov, FVE_threshold)
%=======
%Input:
%=======
% xcov:             ngrid * ngrid matrix of smooth covariance surface
% FVE_threshold:    a positive number that is between 0 and 1
%                   this is the threshold value such that the number of 
%                   selected PC will explain at least this amount of
%                   percentage of total variation
%=======
%Output:
%=======
% no_opt:          positive integer denotes the best number of PC chosen from FVE method
% FVE   :          ngrid * 1 vector cumulative percentage of variation explained by ngrid PC's
%                  it is a cumulative fraction corresponding to d
% lambda     :     ngrid * 1 vector of positive eigenvalues ordered from largest to smallest obtained
%                  from the singular value decomposition of the smooth covariance matrix xcov                    
% Note that FVE and lambda might be smaller if some of the eigenvalues are negative or complex numbers.
function [no_opt, FVE, lambda] = no_FVE(xcov, FVE_threshold, varargin)

     %[u d eigen] = svds(xcov, size(xcov,1), 'L');  %obtain ngrid * 1 eigenvalues
     opts.disp = 0;
     numGrids = size(xcov,1)-2;
     d = eigs(xcov,numGrids,'lm',opts);          %at most ngrid - 2 eigenvalues can be obtained 
                                              %for nonsymmetric and complex problems
     idx = find(imag(d));                    %remove any imaginary eigenvalues
     if isempty(idx) == 0
        %if length(varargin) > 0
        %  fprintf(1,['Warning: ' num2str(length(idx)) ' eigenvalues are complex and are removed!\n']);  
        %end
        %d = d(find(~imag(d)));
        error('Error: %d eigenvalues are complex. The estimated auto-covariance surface is not symmetric!',length(idx));
     end
     idx = find(d <= 0);                       
     if isempty(idx) == 0
        if length(varargin)>0
          fprintf(1,['Warning: ' num2str(length(idx)) ' real eigenvalues are negative or zero and are removed!\n']);
        end  
        d = d(d > 0);                                  %retain only the positive eigenvalues
     end

     %d = svds(xcov, size(xcov,1), 'L');  %obtain ngrid * 1 eigenvalues
     lambda = d';
     FVE = cumsum(lambda)./sum(lambda);
     no_opt = find(FVE >= FVE_threshold, 1, 'first');
end
