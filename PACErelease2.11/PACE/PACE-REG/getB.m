% function [b] = getB(pc_x, yy, y,t_y, isYFun, K_x, K_y,method)
%=======
%Input:
%=======
%  pc_x:  The first K_x FPC scores of X
%    yy:  The returned values from FPCA. See FPCA() for more details.
% y    :  The original response variable, used when the response is scalar        
% t_y  :  The original input observation time for response
% isYFun:  a logical value, isYFun = 1 (Y is a functional response)
%                           isYFun = 0 (Y is a scalar response)
% K_x:   positive integer; number of principal components of predictor x used in regression
%        must be smaller than or equal to the maximum number of principal components give by FPAC
%        default is the number of pc selected by FPCA
% K_y:   positive integer; number of principal components of response y used in regression
%        must be smaller than or equal to the maximum number of principal components give by FPAC
%        default is the number of pc selected by FPCA
% method:  a string taking the following two values:
%            'YPC': fitting the functional linear regression;
%            'YCO': predicting with a functional linear regression. 
%          Refer to the input FIT of FPCreg for more details.
%=======
%Output:
%=======
%  b:   estimated coefficients of beta in the basis represatation 
%       a K_x*K_y matrix if Y is functional; a 1*K_x vector is Y is scalar 

function [b] = getB(pc_x, yy, y,t_y, isYFun, K_x, K_y,method)
        
    if isYFun == 1
           
        if strcmp(method,'YPC') == 1

            pc_y = getVal(yy,'xi_est');
            pc_y = pc_y(:,1:K_y);     % truncated PC for y
            b = zeros(K_x,K_y);
            for j = 1:K_y
                for k = 1:K_x
                    b(k,j) = pc_y(:,j)'*pc_x(:,k)/(pc_x(:,k)'*pc_x(:,k));
                end
            end
            
        elseif strcmp(method,'YCO') == 1
            
            phi_y = getVal(yy,'phicopy');   %corresponding to distinct time points of t_y
            phi_y = phi_y(:,1:K_y);   % truncated eigenfunction
            out1_y = getVal(yy,'out1copy');
            if isempty(getVal(yy, 'sigma'))
                sigma = 0;
            else
                sigma = max(getVal(yy, 'sigmanew'),getVal(yy,'rho_opt'));
            end
            lambda_y = getVal(yy,'lambda');
            lambda_y = lambda_y(1:K_y);
            
            [X W] = getXW(pc_x, phi_y, out1_y, t_y, lambda_y, sigma, 1);
            X = (cell2mat(X))';
            W = blkdiag(W{:});
            tty = cell2mat(t_y);
            Y = cell2mat(y)'-interp1(out1_y,getVal(yy,'mu'),tty,'spline')';
            b = pinv(X'*W*X)*X'*W*Y;
%             b = pinv(X'*X)*X'*Y;
            b = reshape(b,K_y,K_x)';
            
        end
        
    else
        
            pc_y = y-mean(y);
            b = zeros(1,K_x);
            for k = 1:K_x
                b(k) = pc_y*pc_x(:,k)/(pc_x(:,k)'*pc_x(:,k));
            end

    end

