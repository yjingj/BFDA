% function [rej pv Qsim] = Rtest(res, xx, yy, nsim, alpha)
% Rtest    conducts a hypothesis testing of Quasi R-square Q: Q=0, based on 
%          bootstraped sampling; returns decision binary and p-value, by
%          utilizing the output of FPCreg.
%          this bootstrap p-value is based on resampling from the
%          estimated functional principal scores (this method is fast, but 
%          does not reflect the variability in estimating the FPC
%          scores). For a more accurate bootstrap method, based on bootstrapping 
%          by resampling from the original subject data, use FPCdiag.m with 
%          boostrap = 'pre' option.
% 
% Input:
%     res:   the output from FPCreg
%  xx, yy:   The returned values from FPCA. See FPCA() for more details.
%            Let yy=y if y is scalar.
%    nsim:   number of bootstrap sample size; nsim = 500 by default.
%   alpha:   level of hypothesis test, optional;  alpha is 0.05 by default.
% 
% Output:
%     rej:   1 ==> reject the null hypothesis of Q=0;
%            0 ==> do not reject.
%     pv:    P-value of hypothesis testing Q=0.
function [rej pv Qsim] = Rtest(res, xx, yy, nsim, alpha)
    
    isYFun = getVal(res, 'isYFun');
    method = getVal(res, 'method');
    t_x = getVal(xx,'t');
    K_x = getVal(res,'K_x');
    if isYFun == 1
        K_y = getVal(res,'K_y');
    else
        K_y = [];
    end
    if nargin < 4 || isempty(nsim)
        nsim = 500;
    end
    if nargin < 5 || isempty(alpha)
        alpha=.05;
    end
    
    pc_x = getVal(xx,'xi_est');
    pc_x = pc_x(:,1:K_x);
    n = size(pc_x,1);
    if isYFun == 1
        pc_y = getVal(yy,'xi_est');
        pc_y = pc_y(:,1:K_y);
    end
   
    Qsim=zeros(nsim,1);
    rnd = mysample(1:n, n*nsim*2,1);
    rnd = reshape(rnd, n, nsim, 2);
    if isYFun == 1
        for l=1:nsim
            y = getVal(yy, 'y');
            t_y = getVal(yy, 't');
            rndx = rnd(:, l, 1);
            rndy = rnd(:, l, 2);
            y1 = y(rndy);
            t_y1 = t_y(rndy);
            yy1 = setVal(yy, 'xi_est', pc_y(rndy,:));
            b = getB(pc_x(rndx,:),yy1, y1,t_y1,1,K_x,K_y,method);
            fitted_y = predict(pc_x(rndx,:),t_x(rndx),xx,y1,yy1,t_y1,b,isYFun,K_x,K_y,method);
            Qsim(l) = getQ(fitted_y,yy1,y1,t_y1,1);
        end
    else
        for l=1:nsim
            rndx = rnd(:, l, 1);
            rndy = rnd(:, l, 2);
            b = getB(pc_x(rndx,:),[],yy(rndy),[], 0,K_x, K_y,method);
            fitted_y = predict(pc_x(rndx,:),t_x(rndx),xx,yy(rndy),[],[],b,0,K_x,K_y,method);
            Qsim(l) = getQ(fitted_y,[],yy(rndy),[],0);
        end
    end
    
    pv = sum(Qsim>getVal(res,'Q'))/nsim;
    if pv < alpha
        rej = 1;
            if pv == 0
                fprintf(1,['P-value < ' num2str(1/nsim) '. Reject the null hypothesis of no regression relation.\nX and Y are significantly correlated.\n']);
            else
                fprintf(1,['P-value = ' num2str(pv) '. Reject the null hypothesis of no regression relation.\nX and Y are significantly correlated.\n']);
            end
    else
        rej = 0;
        fprintf(1,['P-value = ' num2str(pv) '. Cannot reject the null hypothesis of no regression relation.\nX and Y are not significantly correlated.\n']);
    end
end

