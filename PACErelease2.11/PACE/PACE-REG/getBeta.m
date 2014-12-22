% function [BETA,b] = getBeta(xx, yy, y, t_y, isYFun, K_x, K_y, method)
%=======
%Input:
%=======
% xx, yy:  The returned values from FPCA. See FPCA() for more details.
% y    :  The original response variable.        
% t_y  :  The original input observation time for response.
% isYFun:  a logical value, isYFun = 1 (Y is a functional response)
%                           isYFun = 0 (Y is a scalar response)
% K_x:   positive integer; number of principal components of predictor x used in regression
%        must be smaller than or equal to the maximum number of principal components give by FPAC
%        default is the number of pc selected by FPCA.
% K_y:   positive integer; number of principal components of response y used in regression
%        must be smaller than or equal to the maximum number of principal components give by FPAC
%        default is the number of pc selected by FPCA.
% method:  a string taking the following two values:
%            'YPC': fitting the functional linear regression;
%            'YCO': predicting with a functional linear regression. 
%          Refer to the input FIT of FPCreg for more details.
%=======
%Output:
%=======
% BETA: a struct object that contains:
%       beta: a ngrid_x*ngrid_y matrix; estimated beta function at given grid points  
%       grid_x : 1*ngrid_x vector, it contains ngrid of distinct time points of x function
%                where ngrid is defined in setOptions()
%       grid_y : 1*ngrid_y vector, it contains ngrid of distinct time points of y function
%                where ngrid is defined in setOptions()
%                when y is a scalar, it is []
%       Note: BETA is not calculated when method = 'FCO'.
%  b:   estimated coefficients of beta in the basis represatation 
%       a K_x*K_y matrix if Y is functional; a 1*K_x vector is Y is scalar 

function [BETA,b] = getBeta(xx, yy, y, t_y, isYFun, K_x, K_y, method)

    pc_x = getVal(xx,'xi_est');
    pc_x = pc_x(:,1:K_x);       % truncated PC for x        
    [b] = getB(pc_x, yy, y,t_y, isYFun, K_x, K_y,method);

    phi_x = getVal(xx,'eigencopy');    %corresponding to distinct time points of t_x
    phi_X = phi_x(:,1:K_x);      % truncated eigenfunction

    if strcmp(method,'YPC') == 1 || isYFun == 0
        
        if isYFun == 1
            phi_y = getVal(yy,'eigencopy');   %corresponding to distinct time points of t_y
            phi_Y = phi_y(:,1:K_y);   % truncated eigenfunction
            S = size(phi_X,1);
            T = size(phi_Y,1);
            beta = zeros(S,T);
            for t = 1:T
                for s = 1:S
                    beta(s,t) = phi_X(s,:)*b*phi_Y(t,:)';
                end
            end
        else
            beta = b*phi_X';
        end

        grid_x = getVal(xx,'out21copy');
        if isempty(yy)
            grid_y = [];
        else
            grid_y = getVal(yy,'out21copy');
        end
        BETANames = {'beta','grid_x','grid_y'};
        BETA = {beta,grid_x,grid_y,BETANames};
        
    elseif strcmp(method,'YCO') == 1
        
        BETA = [];

    end
