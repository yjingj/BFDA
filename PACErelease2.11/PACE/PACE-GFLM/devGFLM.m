% calculate Pearson Chi-square and deviance
% D1 - Pearson Chi-sqaure
% D2 - deviance (only for parametric model with canonical link function)

function [D1,D2] = devGFLM(Y,mu,sigma2,vartyp);

D1 = sum((Y-mu).^2./sigma2);
if vartyp==0
    D2=NaN;
elseif vartyp==1
    D2=D1;
elseif vartyp==2
    D2=2*sum(Y.*log((Y+0.000001)./mu)+(1-Y).*log((1-Y+0.000001)./(1-mu)));
elseif vartyp==3
    D2=2*sum(Y.*log(Y./mu)+(mu-Y));
elseif vartyp==4
    D2=2*sum(Y.*(mu-Y)+log(mu./Y));
end
   
end