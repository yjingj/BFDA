% function for updating b, the regression coefficents of FPC scores or X
% (if PACE not used)
% INPUTs: Xreg a n*p matrix of the FPC scores or X (if PACE not used);
% Y,eta,mu,der,sigma as in the inputs of iterGFLM
% OUTPUT: b a 1*p vector

function b = betaupd(Xreg,Y,eta,mu,der,sigma2)

n=size(Xreg,1);
p=size(Xreg,2);
tmptindx=find(sigma2>0);
if length(tmptindx)==0
    tmpt=der.^2./(var(Y)*ones(1,n));
    clear tmptindx
    tmptindx=1:n;
elseif length(tmptindx)==n  
    tmpt=der.^2./sigma2;
else
    tmpt=zeros(1,n);
    tmpt(tmptindx)=der(tmptindx).^2./sigma2(tmptindx);
end
if min(tmpt(tmptindx))==0
    tmptindx1=find(tmpt(tmptindx)==0);
    tmptindx2=find(tmpt(tmptindx)~=0);
    tmpt(tmptindx(tmptindx1))=repmat(min([0.01,tmpt(tmptindx(tmptindx2))]),[1,length(tmptindx1)]);
end
W=repmat(max(tmpt),[1,n]);
W(tmptindx)=tmpt(tmptindx);
W=diag(W);
clear tmptindx tmptindx1 tmptindx2;
tmptindx=find(der==0);
if length(tmptindx)>0
    tmptindx1=find(der~=0);
    der(tmptindx)=(2*(rand(1,length(tmptindx))>0.5)-1)*max(abs(der(tmptindx1)));
end
b=(pinv(Xreg'*W*Xreg)*Xreg'*W*(eta+(Y-mu)./der)')';

end

