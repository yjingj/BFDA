function Dhval = mongrad_ODE(x, Wfd, EPS)
% MONFN_ODE evaluates a monotone function  h(x) = (D^{-1} exp Wfd)(x)
%  where  D^{-1} means taking the indefinite integral.
%  The interval over which the integration takes places is defined in
%       the basis object in WFD.
%  In contrast to the older MONFN, this function approximates the 
%  solution by solving the initial value differential equation
%  corresponding to the indefinite integral defining the function.

%  Arguments:
%  X         ... Argument values at which function and derivatives are 
%                evaluated
%  WFD       ... A functional data object
%  EPS       ... Relative error needed for convergence. Default 1e-8
%
%  Returns:
%  HVAL      ... value of h at input argument array X in first column.

%  Last modified 11 November 2007

if nargin < 2,  error('There are less than two arguments');  end

%  set some constants

if nargin < 3,  EPS = 1e-6;  end

wbasis = getbasis(Wfd);
wrange = getbasisrange(wbasis);
width  = wrange(2) - wrange(1);
nbasis = getnbasis(wbasis);

odeoptions = odeset('RelTol', EPS, 'AbsTol', EPS, ...
                    'MaxStep',     1e-1*width,        ...
                    'InitialStep', 1e-1*width);

[xval, Dhval] = ode23(@Dygradfn, x, zeros(nbasis,1), odeoptions, Wfd);

%  --------------------------------------------------------------------

function Dygradval = Dygradfn(xval, yval, Wfd)
Wval  = eval_fd(xval, Wfd);
wbasis = getbasis(Wfd);
wbasisvec = eval_basis(xval, wbasis)';
Dygradval = wbasisvec.*exp(Wval);
