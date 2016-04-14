function basismat = QWM(u, nderiv)
%QWM computes values of the three basis functions for log-derivative of
%  the modified Weibull quantile function or its derivative.  
%  These are:
%     { 1,  -log(1-u), log(1 - log(1-u)) } = {F1, F2, F3}
%  They are defined over [0, Umax], Umax < 1.
%  Argument U is checked for satisfying this condition.

%  Arguments are as follows:
%  U       ... A vector of values at which the basis functions are to be
%              evaluated
%  NDERIV  ... highest order derivative.  Cannot exceed 3
%  Return is a matrix with length(U) rows and 3 columns containing
%  the values of the QMW basis functions or their derivatives.

%  last modified 11 June 2009

if nargin < 2, nderiv = 0;  end

if nderiv > 3
    error('NDERIV is greater than 3.');
end

if max(u) >= 1
    error('Maximum argument value greater than or equal to 1.');
end

if min(u) < 0
    error('Minimum argument value less than 0.');
end

n  = length(u);
F1 = ones(n,1);
F2 = -log(1-u);
F3 = log(1 + F2);
basismat = zeros(n,3);

if nderiv == 0
    basismat(:,1) = F1;
    basismat(:,2) = F2;
    basismat(:,3) = F3;
end
 
if nderiv == 1
    basismat(:,2) = 1./(1-u);
    basismat(:,3) = 1./((1-u).*(1+F2));
end

if nderiv == 2
    basismat(:,2) = 1./(1-u).^2;
    basismat(:,3) = F2./((1-u).^2.*(1+F2).^2);
end

if nderiv == 3
    basismat(:,2) = 2./(1-u).^3;
    basismat(:,3) = (2.*F2.^2 + F2 + 1)./((1-u).^3.*(1+F2).^3);
end



