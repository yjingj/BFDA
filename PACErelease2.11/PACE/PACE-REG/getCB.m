% function [cb grid_x grid_y] = getCB(res, xx, yy, K_x, K_y, isYFun, alpha)
% getCB    returns 100(1-alpha)% confidence bands of regression function
%          Beta and plot, by utilizing the output of FPCreg.
% Input:
%       res: the output from FPCreg
% xx, yy:  The returned values from FPCA. See FPCA() for more details.
%           Let yy=y if y is scalar.
%      K_x: number of principal components of predictor X used in regression
%      K_y: number of principal components of response Y used in regression
%   alpha:  level of confidence bands, optional; 100(1-alpha)% confidence
%           bands will be created.  alpha is .05 by default.
% Output:
%       cb :     a ngrid_x*ngrid_y*2 matrix; 
%                cb(:,:,1) is the estimated lower confidence band at given grid points;
%                cb(:,:,2) is the estimated upper confidence band at given grid points;
%       grid_x : 1*ngrid_x vector, it contains ngrid of distinct time points of x function
%                where ngrid is defined in setOptions()
%       grid_y : 1*ngrid_y vector, it contains ngrid of distinct time points of y function
%                where ngrid is defined in setOptions()
%                when y is a scalar, it is []

function [cb grid_x grid_y] = getCB(res, xx, yy, K_x, K_y, alpha)
    
    cb = []; grid_x = []; grid_y = [];
    isYFun = getVal(res, 'isYFun');
    method = getVal(res, 'method');
    if strcmp(method,'YPC') == 1
        nsim = 1000;
    elseif isYFun == 1
%         nsim = 200;
        fprintf(1,'Warning: Diagnostics of the fitted functional linear regression model is available only for the case FIT = 0.\n');
        return;
    else
        nsim = 1000;
    end
    if isempty(K_x) 
        K_x = getVal(res,'K_x');
    end
    if isempty(K_y)&&(isYFun==1) 
        K_y = getVal(res,'K_y');
    end
    if nargin<7||isempty(alpha)
        alpha=0.05;
    end
    if alpha <=0 || alpha >= 1
       fprintf(1,'Warning: alpha must be between 0 and 1! Reset alpha = 0.05 now!\n');
    end
    
    grid_x = getVal(xx,'out21copy');    
    pc_x = getVal(xx,'xi_est');
    pc_x = pc_x(:,1:K_x);
    phi_x = getVal(xx,'eigencopy');    %corresponding to distinct time points of t_x 
    phi_x = phi_x(:,1:K_x);      % truncated eigenfunction
    S = length(phi_x(:,1));
    n = length(pc_x(:,1));
    grid_y = [];
    if isYFun == 1
        grid_y = getVal(yy,'out21copy');
        pc_y = getVal(yy,'xi_est');
        pc_y = pc_y(:,1:K_y);
        phi_y = getVal(yy,'eigencopy');   %corresponding to distinct time points of t_y
        phi_y = phi_y(:,1:K_y);   % truncated eigenfunction
        T = length(phi_y(:,1));
    end
    
    rnd = mysample(1:n, n*nsim,1);
    rnd = reshape(rnd, n, nsim);
    figure;
    if isYFun == 1
        y = getVal(yy, 'y');
        t_y = getVal(yy, 't');
        beta2=zeros(S,T,nsim);
        for l=1:nsim
            rnd1 = rnd(:,l);
            yy1 = setVal(yy, 'xi_est', pc_y(rnd1,:));
            b = getB(pc_x(rnd1,:), yy1,  y(rnd1), t_y(rnd1), isYFun, K_x, K_y, method);
            for t = 1:T
                for s = 1:S
                    beta2(s,t,l) = phi_x(s,:)*b*phi_y(t,:)';
                end
            end
        end
        cb=zeros(length(grid_x),length(grid_y),2);
        for j=1:S
            for k=1:T
                cb(j,k,:)=quantile(beta2(j,k,:),[alpha/200 1-alpha/200]);
            end
        end
        
        surf(grid_x,grid_y,cb(:,:,1));
        hold on;
        surf(grid_x,grid_y,cb(:,:,2));
        surf(grid_x,grid_y,getVal(getVal(res,'BETA'),'beta'));
        hold off;
        xlabel('x');
        ylabel('y');
        title('Confidence bands for Beta');
    else
        beta1=zeros(S,nsim);
        yy = yy - mean(yy);
        for l=1:nsim
            rndx = rnd(:,l);
            b = getB(pc_x(rndx,:), [], yy(rndx), [], isYFun, K_x, 1, method);
            beta1(:,l) = b*phi_x';
        end
        cb =zeros(S,2);
        for j=1:S
            cb(j,:)=quantile(beta1(j,:),[alpha/2 1-alpha/2]);
        end
        
        plot(grid_x, cb(:,1),'--');
        hold on;
        plot(grid_x, cb(:,2),'--');
        plot(grid_x, getVal(getVal(res,'BETA'),'beta'));
        hold off;
        xlabel('x');
        title('Confidence bands for Beta');
    end
end
