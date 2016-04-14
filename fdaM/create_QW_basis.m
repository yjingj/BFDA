function basisobj = create_QW_basis(rangeval)
%  CREATE_QWM_BASIS  Creates a basis for the log-derivative W of the
%  Weibull quantile function.  
%  These basis functions are:
%     { 1,  -log(1-u), log(-log(1-u)) }
%  They are defined over [Umin, Umax], Umin > 0,  Umax < 1.
%  Argument RANGEVAL is checked for satisfying this condition.
%  RANGEVAL ... Range over which the quantile function is defined.
%               The upper limit MUST be less than 1, and the lower
%               limit MUST be greater than 0.  
%  Return:
%  BASIS.FD  ... a functional data basis object of type 'QW'
%

%  Last modified 11 June 2009

%  Check rangeval

if length(rangeval) == 1
    rangeval = [0, rangeval];
end

if rangeval(1) <= 0
    error('RANGEVAL(1) <= 0.');
end

if rangeval(2) >= 1
    error('RANGEVAL(2) >= 1');
end

type        = 'QW';
params      = [];
nbasis      = 3;
dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

