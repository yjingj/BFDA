function [hik,f]=rthik(t,curvei,curvek,nknots,lambda)
% find the pairwise warping function between curve i and k, curve k is the reference curve
% input
% t,curvei,curvek all 1*m vector, give the time points, vectors for curve i and curve j
% nknots, number of knots used when estimating pairwise warping functions
% lambda, penalty parameter used when estimating pairwise warping functions
% output
% hik, the estimated pairwise warping function
% f, value of the minimizer function at solution hik

warning off all 
temp=linspace(min(t),max(t),nknots+2);
struct=temp(2:nknots+1); 
lb=zeros(size(struct))+0.0001';
ub=max(t)*ones(size(struct'))-0.0001;
b=zeros(nknots-1,1)-0.0002;
C=[eye(nknots-1),zeros(nknots-1,1)]-[zeros(nknots-1,1),eye(nknots-1)];
initial=struct;
[htemp,f]=fmincon(@rttemp1,initial,C,b,[],[],lb,ub,[],optimset('Diagnostics','off','display','off'));

    function f=rttemp1(htemp)
        hik=fnval(spapi(2,[min(t),struct,max(t)],[min(t),htemp,max(t)]),t);
        f=sum((fnval(curvei,hik)-fnval(curvek,t)).^2);
        penal=hik-t;
        f=f+lambda*sum(penal.^2);
    end

end

