% This function is to compute the leave-one-point-out prediction error, 
% mainly for choosing the overdispersion parameter sigma (usually in conjunction 
% with parametric variance function).
% The LOPO is w.r.t. removing Y_{ij}, while the model components, such as mean 
% and eigenfunctions are unchanged.
function [ypred,indpred,yprederr]=grm_predict_lopo(yo,to,mux,muy,lambdax,phix,sigma,...
    vf,family,param,ind_pred,out1)

 % if sigma==1, the variance function v(.) is usually nonparametrically estimated and given by vf(2*m); 
% if the variance function is of a known parametric form, then it is evaluated by varf, and vf==[]; 

ntemp = [];
for i = 1:length(ind_pred)
    ntemp(i) = length(to{ind_pred(i)});
end
ind_pred1 = ind_pred(ntemp>1);
ypred = {}; yprederr = [];

for i = 1:length(ind_pred1)
    get = ind_pred1(i);    
    muytemp = interp1(out1,muy,to{get},'spline');
    muxtemp = interp1(out1,mux,to{get},'spline');
    phii = interp1(out1',phix,to{get}','spline');
    if isempty(vf)
        error0 = sigma*diag(varf(muytemp,family,param));
    else
        temp2 = glink(muxtemp,family,param);
        error0 = diag(interp1([0 vf(1,:) 1],[0 vf(2,:) 0],temp2,'linear'));
    end
    
    for j = 1:length(to{get})
        indt = 1:length(to{get});
        indt(j) = [];
        phiit = phii(indt,:);
        xi_predt = (diag(lambdax)*phiit'*diag(glink_der(muxtemp(indt),family,param))*pinv(diag(glink_der(muxtemp(indt),family,param))*phiit*...
                 diag(lambdax)*phiit'*diag(glink_der(muxtemp(indt),family,param))+error0(indt,indt))*(yo{get}(indt)-muytemp(indt))')';
        temp1 = phii(j,:);
        temp2 = muxtemp(j)+xi_predt*temp1';
        ypred{get}(j) = glink(temp2,family,param);             
    end         
end

for i = 1:length(ind_pred1)
    get = ind_pred1(i);
    if isempty(vf)
        yprederr(i) = sum((yo{get}-ypred{get}).^2./(varf(ypred{get},family,param)));
        %yprederr(i)=sum((yo{get}-ypred0{get}).^2);        
    else
        temp3 = interp1([0 vf(1,:) 1],[0 vf(2,:) 0],ypred{get},'linear');
        yprederr(i) = sum((yo{get}-ypred{get}).^2./temp3);
        %yprederr(i)=sum((yo{get}-ypred0{get}).^2); 
    end
end
indpred = ind_pred1;
