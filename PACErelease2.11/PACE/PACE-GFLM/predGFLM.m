%===========
%Description:
%===========
%
%            This program is needed to obtain predicted values at new
%            covariate levels. Estimated link function and regression
%            coefficient function must be given on a dense grid. Linear
%            interpolation is used to obtain the values of these functions
%            as needed to obtain the predicted values.
%  
%======
%Usage:
%======
%
% function [newY,neweta] = predGFLM(newX,newT,beta,tgrid,alpha,ggrid,gfct)
%
%======
%Input: 
%======
%      newX:      n1*m matrix or 1*n1 cell array storing the new predictors 
%      newT:      1*m vector or 1*n1 cell array giving the time points 
%      beta:      1*m vector giving the regression parameter function 
%      tgrid:     1*m vector givng the time grid for beta
%      alpha:     estimated intercept (0 if link function is unknown) 
%      ggrid:     1*100 vector giving the time grid for gfct from iterGFLM
%      gfct:      1*100 vector giving the evaluated link function on time
%                 grid ggrid from iterGFLM
%
%=======
%Output:  
%=======  
%      newY:      1*n1 vector giving the predicted g(neweta)
%      neweta:    1*n1 vector giving the predicted f(E(Y|mu)) for newX   

function [newY,neweta] = predGLM(newX,newT,beta,tgrid,alpha,ggrid,gfct)

if ~iscell(newT)&all(newT-tgrid==0)
    neweta=trapz(newT',newX'.*(beta'*ones(1,size(newX,1))))+alpha;
elseif ~iscell(newT)
    bound1=max(min(newT{k}),min(tgrid));
    bound2=min(max(newT{k}),max(tgrid));
    tgridnew=bound1:(bound2-bound1)/99:bound2;
    for k=1:size(newX,1)
        neweta(k)=trapz(tgridnew,fnval(spapi(2,newT,newX(k,:)),tgridnew).*fnval(spapi(2,tgrid,beta),tgridnew))+alpha;
    end
else
    for k=1:size(newX,1)
        tgridk=[];
        bound1=max(min(newT{k}),min(tgrid));
        bound2=min(max(newT{k}),max(tgrid));
        tgridk=bound1:(bound2-bound1)/99:bound2;
        neweta(k)=trapz(tgridk,fnval(spapi(2,newT{k},newX{k}),tgridk).*fnval(spapi(2,tgrid,beta),tgridk))+alpha;
    end
end
newY=fnval(spapi(2,ggrid,gfct),neweta);

end

