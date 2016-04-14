function basisobj = create_QWM_basis(rangeval)
%  CREATE_QWM_BASIS  Creates a basis for the log-derivative W of the
%  modified Weibull quantile function.  
%  These basis functions are:
%     { 1,  -log(1-u), log(1 - log(1-u)) }
%  They are defined over [0, Umax], Umax < 1.
%  Argument RANGEVAL is checked for satisfying this condition.

%  Arguments:  None required
%  RANGEVAL ... Range over which the quantile function is defined.
%               The upper limit MUST be less than 1, and the lower
%               limit will normally be 0.  If RANGEVAL is a single
%               number, this will be taken as the upper limit, and the
%               lower limit as 0.
%  Return:
%  BASIS.FD  ... a functional data basis object of type 'QWM'
%

%  Last modified 11 June 2009

%  Check rangeval

if length(rangeval) == 1
    rangeval = [0, rangeval];
end

if rangeval(1) < 0
    error('RANGEVAL(1) < 0.');
end

if rangeval(2) >= 1
    error('RANGEVAL(2) >= 1');
end

type        = 'QWM';
params      = [];
nbasis      = 3;
dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

