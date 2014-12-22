%function [beta1] = getDynamic(yy)
%This function calculate the dynamic transfer function,
%\beta_1(t) for nder = 1.
%
%=====
%Input
%=====
%
%  yy : returned object from FPCder().
%
%======
%Output
%======
%  beta: 1*N vector, estimated dynamic transfer function evaluated 
%        at out1 or [] if none of the estimated derivative is related
%        to nder = 1.
%
%example:
%yy = FPCder(y,t,p);
%beta1 = getDynamic(yy)

function [beta1] = getDynamic(yy)

  ops = getVal(yy,'ops');
  nder = getVal(ops,'nder');

  if isempty(find(nder ==1)) == 0

    lambda = getVal(yy,'lambda');
    phi = getVal(yy,'phi');
    phi0 = phi{find(nder == 0)};
    phi1 = phi{find(nder == 1)};
    xcov0 = phi0*diag(lambda)*phi0';
    bottom = diag(xcov0);
    xcov1 = phi0*diag(lambda)*phi1';
    top = diag(xcov1);
    beta1 = top./bottom;
    if size(beta1,2) == 1
      beta1 = beta1';
    end
  else
    fprintf(1, 'Warning: only nder = 1 is supported for beta function calculation!\n');
    beta1 = [];
  end

end
