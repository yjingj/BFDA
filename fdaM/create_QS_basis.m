function basisobj = create_QS_basis(rangeval, nbasis)
%  CREATE_QS_BASIS  Creates a basis for the log-derivative W of the
%  a quantile function that is a function of surprisal s.  
%  These basis functions are:
%     { 1,  s, log(1 + s), log[1 + log(1 + s)], ... }
%  They are defined over [0, Smax], but Smax is not actually used.
%  Argument RANGEVAL is checked for satisfying this condition.

%  Arguments:  None required
%  RANGEVAL ... Range over which the quantile function is defined.
%               The upper limit MUST be less than 1, and the lower
%               limit will normally be 0.  If RANGEVAL is a single
%               number, this will be taken as the upper limit, and the
%               lower limit as 0.
%  Return:
%  BASIS.FD  ... a functional data basis object of type 'QS'
%

%  Last modified 9 December 2009

%  set default NBASIS

if nargin < 2,  nbasis = 2;  end

%  Check rangeval

if length(rangeval) == 1
    rangeval = [0, rangeval];
end

if rangeval(1) < 0
    error('RANGEVAL(1) < 0.');
end

type        = 'QS';
params      = [];
dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

