function basismat = QS(s, nbasis, nderiv)
%QS computes values of the basis functions for log-derivative of
%  the quantile function as a function of surprisal s
%  or its derivative.  
%  These are:
%     { 1,  s, log(1+s), , log[1+log(1+s)], ... } 
%  They are defined over [0, Smax], but Smax is never actually
%  used
%  Argument S is checked for being nonnegative

%  Arguments are as follows:
%  S       ... A vector of values at which the basis functions are to be
%              evaluated
%  NBASIS  ... Number of basis functions required
%  NDERIV  ... highest order derivative.  
%  Return is a matrix with length(U) rows and NBASIS columns containing
%  the values of the QS basis functions or their derivatives.

%  last modified 9 December 2009

if nargin < 3, nderiv = 0;  end

if nderiv > 3
    error('Derivatives up to order 3 only are allowed.');
end

if nargin < 2, nbasis = 2;  end

if min(s) < 0
    error('Minimum argument value less than 0.');
end

n  = length(s);

%  compute basis functions 

basismat = zeros(n,nbasis);

basismat(:,1) = ones(n,1);
if nbasis > 1
    basismat(:,2) = s;
end
if nbasis > 2
    for ibasis=3:nbasis
        basismat(:,ibasis) = log(1 + basismat(:,ibasis-1));
    end
end

%  compute derivative of basis functions

if nderiv >= 1
    D1basismat = zeros(n,nbasis);
    if nbasis > 1, D1basismat(:,2) = ones(n,1);  end
    if nbasis > 2
        for ibasis = 3:nbasis
            D1basismat(:,ibasis) = D1basismat(:,ibasis-1)./ ...
                (1 + basismat(:,ibasis-1));
        end
    end
    if nderiv == 1
        basismat = D1basismat;
        return;
    end
end

if nderiv >= 2
    D2basismat = zeros(n,nbasis);
    if nbasis > 2
        D2basismat(:,3) = -1./(1+s).^2;
    end
    if nbasis > 3
        for ibasis = 4:nbasis
            D2basismat(:,ibasis) = - D1basismat(:,ibasis-1)./      ...
                                     (1 + basismat(:,ibasis-1)).^2 ...
                                   + D2basismat(:,ibasis-1)./      ...
                                     (1 + basismat(:,ibasis-1));
        end
    end
    if nderiv == 2
        basismat = D2basismat;
        return;
    end
end

if nderiv == 3
    D3basismat = zeros(n,nbasis);
    if nbasis > 2
        D3basismat(:,3) = 2./(1+s).^3;
    end
    if nbasis > 3
        for ibasis = 4:nbasis
            D3basismat(:,ibasis) =  2.*D1basismat(:,ibasis-1).^2./  ...
                                     (1 + basismat(:,ibasis-1)).^3 ...
                                   -2.*D2basismat(:,ibasis-1)./      ...
                                     (1 + basismat(:,ibasis-1)).^2 ...
                                      +D3basismat(:,ibasis-1)./ ...
                                     (1 + basismat(:,ibasis-1));
        end
    end
    basismat = D3basismat;
end

