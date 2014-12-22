% Choose number of eigenfunctions K using AIC/BIC based on quasi-likelihood when 
% parametric variance function and overdispersion parameter are used. 
% This is usually used in iterative selection procedure in conjunction with 
% the function overdisp.m to choose K and sigma iteratively.
function [no_opt,cret_val] = noeig_param(cret,y,t,mux,muy,lambdax, phix,sigma,family,param,out1)

    ncohort = length(y); 
    ind_pred = 1:ncohort; 
    N = length(cell2mat(y)); 

    maxk = length(lambdax);
    ql = zeros(1,maxk);
    for k = 1:maxk
        [xi_pred,xpred,ypred,xpred0,ypred0,yprederr] = grm_predict(y,t,mux,muy,lambdax(1:k),phix(:,1:k),sigma,...
                               [],family,param,ind_pred,out1,[]);
        for m = 1:ncohort
            ql(k) = ql(k)+sum(QL_param(ypred0{m},y{m},family,param,sigma));
        end
    %     aic(k)=-2*ql(k)+2*k;
    %     bic(k)=-2*ql(k)+log(N)*k;
    end

    if strcmp(cret,'AIC1') == 1
        cret_val = -2*ql+2*(1:maxk);
    elseif strcmp(cret,'BIC1') == 1
        cret_val = -2*ql+log(N)*(1:maxk);
    end
    no_opt = find(cret_val==min(cret_val));


