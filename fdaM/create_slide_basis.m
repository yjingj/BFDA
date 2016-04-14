function basisobj = create_slide_basis(rangeval, breaks, ratevec)
%CREATE_SLIDE_BASIS Creates a consisting of step functions multiplied
%  exponential decays over the intervals between steps.
%  Arguments ...
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  BREAKS   ... also called knots, these are a strictly increasing sequence
%               of junction points between piecewise polynomial segments.
%               They must satisfy BREAKS(1) = RANGEVAL(1) and
%               BREAKS(NBREAKS) = RANGEVAL(2), where NBREAKS is the total
%               number of BREAKS.
%  RATEVEC  ... Rate of decay over each interval.  Must be of length
%               length(BREAKS) - 1.  Default values are one.
%  Returns
%  BASISOBJ  ... a functional data basis object

%  last modified 14 March 2010

if nargin > 3
    error('More than two arguments for slide basis.');
end

%  Default basis for missing arguments

if nargin==0
    type     = 'slide';
    rangeval = [0,1];
    params   = 1;
    nbasis   = 1;

    basisobj = basis(type, rangeval, nbasis, params);
    return
end

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

%  set some default values

if nargin < 2 
    params = 1; 
    nbasis = 1;
end

%  Special argument configurations taken care of.  
%  Now go ahead and set up the basis

breaks  = breaks(:);
nbreaks = length(breaks);
nbasis  = nbreaks - 1;

%  check the compatibility of NBASIS, NBREAKS and RANGEVAL

if (nbreaks < 2)
    error ('Number of values in BREAKS less than 2.');
end

%  set default rate values 

if nargin < 3,  ratevec = ones(nbreaks-1,1);  end

ratevec = ratevec(:);

%  The PARAMS field contains only the interior knots; drop end breaks

if nbreaks > 2
    params   = [breaks(2:(nbreaks-1)); ratevec]';
else
    params = [];
end

%  construct basis object

type = 'slide';

basisobj = basis(type, rangeval, nbasis, params);
