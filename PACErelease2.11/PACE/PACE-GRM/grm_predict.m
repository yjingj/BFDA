% This function is to predict the smoothed individual trajectories for 
% the observed non-Gaussian process Y and the latent Gaussian process X, 
% evaluated at output grid out1.
function [xi_pred,xpred,ypred,xpred0,ypred0,yprederr]=grm_predict(yo,to,mux,muy,lambdax,phix,sigma,...
    vf,family,param,ind_pred,out1,tout)

% if sigma==1, the variance function vf(.) is usually nonparametrically estimated and given by vf(2*m); 
% if the variance function is of a known parametric form, then it is evaluated by varf, and vf==[]; 

    xi_pred = [];
    xi_var = {};
    
    for i = 1:length(ind_pred)
        
        get = ind_pred(i);
        muytemp = interp1(out1,muy,to{get},'spline');
        muxtemp = interp1(out1,mux,to{get},'spline');
        phii = interp1(out1',phix,to{get}','spline');        
        if isempty(vf)
            error0 = sigma*diag(varf(muytemp,family,param));
        else
            temp2 = glink(muxtemp,family,param);
            error0 = diag(interp1([0 vf(1,:) 1],[0 vf(2,:) 0],temp2,'linear'));
        end
        xi_pred(i,:) = (diag(lambdax)*phii'*diag(glink_der(muxtemp,family,param))*pinv(diag(glink_der(muxtemp,family,param))*phii*...
            diag(lambdax)*phii'*diag(glink_der(muxtemp,family,param))+error0)*(yo{get}-muytemp)')';
        xi_var{get} = diag(lambdax)-diag(lambdax)*phii'*diag(glink_der(muxtemp,family,param))*...
            pinv(diag(glink_der(muxtemp,family,param))*phii*diag(lambdax)*phii'*diag(glink_der(muxtemp,family,param))+error0)...
            *(diag(lambdax)*phii'*diag(glink_der(muxtemp,family,param)))';
    end

    xpred = []; xpred0 = {}; 
    ypred = {}; ypred0 = {}; yprederr = [];
    
    if ~isempty(tout)
        muxtemp=interp1(out1,mux,tout,'spline');
        xpred = repmat(muxtemp,length(ind_pred),1)+xi_pred*phix';
        xpred = num2cell(xpred,2)';
        for i = 1:length(ind_pred)
            ypred{i} = glink(xpred{i},family,param);
            %var0=diag(phii*xi_var{get}*phii');
            %rand_int{get}=[xi_est(i,:)*phii'-1.96*sqrt(var0');...
            %               xi_est(i,:)*phii'+1.96*sqrt(var0')];
            %ypred_int{get}=[mutemp+rand_int{get}(1,:);mutemp+rand_int{get}(2,:)];
            %rand_intsm{get}=[xi_est(i,:)*phii'-sqrt(chi2inv(0.95,no_opt))*sqrt(var0');...
            %                 xi_est(i,:)*phii'+sqrt(chi2inv(0.95,no_opt))*sqrt(var0')];
            %ypred_intsm{get}=[mutemp+rand_intsm{get}(1,:);mutemp+rand_intsm{get}(2,:)];
        end
    end

    for i = 1:length(ind_pred)
        get = ind_pred(i);
        muxtemp=interp1(out1,mux,to{get},'spline');
        phii = interp1(out1',phix,to{get}','spline');
        xpred0{get} = muxtemp+xi_pred(i,:)*phii';
        ypred0{get} = glink(xpred0{get},family,param);
        if isempty(vf)
            yprederr(i)=sum((yo{get}-ypred0{get}).^2./(varf(ypred0{get},family,param)));
            %yprederr(i)=sum((yo{get}-ypred0{get}).^2);
        else
            temp3=interp1([0 vf(1,:) 1],[0 vf(2,:) 0],ypred0{get},'linear');
            yprederr(i)=sum((yo{get}-ypred0{get}).^2./temp3);
        end
    end

    
