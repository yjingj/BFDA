function monommat = monomial(evalarg, exponents, nderiv, argtrans)
%  MONOMIAL Values of monomials, or their derivatives.
%  The powers of EVALARG are the NBASIS nonnegative integers in EXPONENTS.
%  The default is 1, meaning EVALARG itself.
%  Arguments are as follows:
%  EVALARG   ... array of values at which the polynomials are to
%                evaluated
%  EXPONENTS ... the set of powers of EVALARG to be used
%  NDERIV    ... order of derivative to be returned.
%  Return is:
%  A matrix with length(EVALARG) rows and NBASIS columns containing
%    the values of the monomials or their derivatives
%  Warning!  There is no protection here against the consequences of using
%     large numbers in evalarg.  If monomials need normalizing, this 
%     must be done externally.

%  last modified 7 July 2014

% set default arguments

if nargin < 4
    argtrans = [0,1];
end

if nargin < 3
    nderiv = 0;
end

if nargin < 2
    exponents = 1;
end

evalarg = evalarg(:);
evalarg = (evalarg - argtrans(1))./argtrans(2);
n = length(evalarg);

nbasis = length(exponents);

%  check whether exponents are nonnegative integers

for ibasis=1:nbasis
    if exponents(ibasis) - round(exponents(ibasis)) ~= 0
        error('An exponent is not an integer.');
    end
    if exponents(ibasis) < 0
        error('An exponent is negative.');
    end
end

% check if there are duplicate exponents

if min(diff(sort(exponents))) == 0
    error('There are duplicate exponents.');
end

monommat = zeros(n,nbasis);

if nderiv == 0
    %  use the recursion formula to compute monomnomial values
    for ibasis=1:nbasis, monommat(:,ibasis) = evalarg.^exponents(ibasis); end
else
    for ibasis=1:nbasis
        degree = exponents(ibasis);
        if nderiv <= degree 
            fac = degree;
            for ideriv=2:nderiv
                fac = fac*(degree-ideriv+1);
            end
            monommat(:,ibasis) = fac.*evalarg.^(degree-nderiv);
        end
    end
end

