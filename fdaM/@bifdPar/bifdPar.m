function bifdParobj = bifdPar(bifdobj, Lfdobjs, Lfdobjt, ...
                              lambdas, lambdat, estimate, penmats, penmatt)
% Sets up a functional parameter object.
%  Arguments:
%  BIFDOBJ  ... A bivariate functional data object.  The basis for this  
%               object defines the bivariate functional parameter.
%               When an initial value is required for iterative 
%               estimation of a bivariate functional parameter, the 
%               coefficients will give initial values for the iteration.
%  LFDOBJS  ... A linear differential operator value or a derivative
%               value for penalizing the roughness of the object
%               with respect to the first argument s.
%               By default, this is 2.
%  LFDOBJT  ... A linear differential operator value or a derivative
%               value for penalizing the roughness of the object
%               with respect to the second argument t.
%               By default, this is 2.
%  LAMBDAS  ... The penalty parameter controlling the smoothness of the 
%               estimated parameter with respect to the first  argument s.  
%               By default this is 0.
%  LAMBDAT  ... The penalty parameter controlling the smoothness of the 
%               estimated parameter with respect to the second argument t.  
%               By default this is 0.
%  ESTIMATE ... If nonzero, the parameter is estimated; if zero, the
%               parameter is held fixed at this value.
%               By default, this is 1.
%  PENMATS  ... A roughness penalty matrix for variation over s
%  PENMATT  ... A roughness penalty matrix for variation over t

%  Last modified 23 November 2012

superiorto('double', 'struct', 'cell', 'char', 'inline', 'basis');

%  check BIFDOBJ

if strcmp(class(bifdobj), 'bibasis')
    [sbasis, tbasis]  = getbasis(bifdobj);
    snbasis = getnbasis(sbasis);
    tnbasis = getnbasis(tbasis);
    bifdobj = bifd(zeros(snbasis,tnbasis),sbasis,tbasis);
end

if ~strcmp(class(bifdobj), 'bifd')
    error('BIFDOBJ is not a bivariate functional data object.');
end

%  set some default argument values

if nargin < 8;  penmatt   = [];          end
if nargin < 7;  penmats   = [];          end
if nargin < 6;  estimate  = 1;           end
if nargin < 5;  lambdat   = 0;           end
if nargin < 4;  lambdas   = 0;           end
if nargin < 3;  Lfdobjt   = int2Lfd(2);  end
if nargin < 2;  Lfdobjs   = int2Lfd(2);  end

%  check the linear differential operators

Lfdobjs = int2Lfd(Lfdobjs);
Lfdobjt = int2Lfd(Lfdobjt);

if ~isa_Lfd(Lfdobjs)
    error('LFDOBJS is not a linear differential operator object.');
end
if ~isa_Lfd(Lfdobjt)
    error('LFDOBJT is not a linear differential operator object.');
end

%  check the roughness penalty parameters

if ~isnumeric(lambdas)
    error('LAMBDAS is not numeric.');
end
if lambdas < 0
    error('LAMBDAS is negative.');
end
if ~isnumeric(lambdat)
    error('LAMBDAT is not numeric.');
end
if lambdat < 0
    error('LAMBDAT is negative.');
end

if ~isnumeric(estimate)
    error('ESTIMATE is not numeric.');
end

%  set up the bifdPar object

bifdParobj.bifd     = bifdobj;
bifdParobj.Lfds     = Lfdobjs;
bifdParobj.Lfdt     = Lfdobjt;
bifdParobj.lambdas  = lambdas;
bifdParobj.lambdat  = lambdat;
bifdParobj.estimate = estimate;
bifdParobj.penmats  = penmats;
bifdParobj.penmatt  = penmatt;

bifdParobj = class(bifdParobj, 'bifdPar');
