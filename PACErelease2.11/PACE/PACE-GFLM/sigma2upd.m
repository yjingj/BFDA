% function for updating sigma2 in the iteration
% INPUTs: Y,theta,mu,vartyp,kernel as in the inputs of iterGFLM; bw gives
% the bandwidth for smoothing variance function
% OUTPUTs: sigma2 as in the input of iterGFLM; optbw3 gives the GCV
% bandwidth (if not specified) 

function [sigma2,optbw3] = sigma2upd(Y,theta,mu,vartyp,kernel,bw)

n=length(Y);
if vartyp==0
    if length(bw)==0
        optbw3=gcv1d((Y-mu).^2,mu,kernel,[],1,0,2,'off',0,ones(1,n));
    else
        optbw3=bw;
    end
    [tmpt,sigma2]=locpoly(optbw3,kernel,[],1,0,mu,(Y-mu)'.^2,ones(1,n),mu);
    clear tmpt; 
else
    sigma2=varfun(mu,Y,vartyp,theta);
    optbw3=NaN;
end
% if min(sigma2)<=0
%     sigma2=max(sigma2,zeros(1,n)+min(sigma2(sigma2>0)));
% end   
if min(sigma2)<=median(sigma2)*0.01
    sigma2=max(sigma2,zeros(1,n)+min(sigma2(sigma2>median(sigma2)*0.01)));
end
    
end

