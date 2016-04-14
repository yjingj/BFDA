function basisobj = create_monomial_basis(rangeval, nbasis, exponents, ...
                                          argtrans, dropind)
%  CREATE_MONOMIAL_BASIS  Creates a monomial basis:, 
%     phi_1(t) = t^{exponents(1)}, phi_2(t) = t^{exponents(2)}, ...
%  The exponents in this version must be nonnegative integers
%  Argument:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS    ... number of basis functions
%  EXPONENTS ... Set of powers of EVALARG to be applied.
%  ARGTRANS  ... A vector of two constants for shifting and scaling the 
%                argument range.  In functions MONOMIAL and MONOM, 
%                the EVALARG is transformed to 
%                [evalarg-argtrans(1))./evalarg(2);
%                Defaults to [0,1];
%  DROPIND  ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Return:
%  BASIS  ... a functional data basis object of type 'monom'

%  last modified 7 July 2014

%  default RANGEVAL

if nargin < 1, rangeval = [0,1];  end

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

%  set up default arguments

if nargin < 2, nbasis = 2;  end
if nargin < 3 || isempty(exponents) 
    exponents = 0:(nbasis-1); 
end
if nargin < 4 || isempty(argtrans) 
    argtrans = [0,1];         
end
if nargin < 5, dropind = [];  end

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

if min(diff(sort(exponents))) <= 0
    error('There are duplicate exponents.');
end

%  check DROPIND

if length(dropind) > 0
    if length(dropind) >= nbasis
        error('Too many index values in DROPIND.');
    end
    dropind = sort(dropind);
    if length(dropind) > 1
        if any(diff(dropind)) == 0
            error('Multiple index values in DROPIND.');
        end
    end
    for i=1:length(dropind);
        if dropind(i) < 1 || dropind(i) > nbasis
            error('An index value is out of range.');
        end
    end
end

parStruct.exponents = exponents;
parStruct.argtrans  = argtrans;

%  set up basis object

type        = 'monom';
params      = parStruct;
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

