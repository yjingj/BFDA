function [fdobj, beta, df, gcv, SSE, dev, ...
          penmat, y2cMap, argvals, y] = ...
                       smooth_GLM(argvals, y, fdParobj, varargin)
%SMOOTH_GLM  Smooths discrete curve represented by basis function
%  expansions fit by penalized least squares.
%
%  Required arguments for this function are:
%
%  ARGVALS  ... A set of argument values, set by default to equally spaced
%               on the unit interval (0,1).
%  Y        ... an array containing values of curves
%               If the array is a matrix, rows must correspond to argument
%               values and columns to replications, and it will be assumed
%               that there is only one variable per observation.
%               If Y is a three-dimensional array, the first dimension
%               corresponds to argument values, the second to replications,
%               and the third to variables within replications.
%               If Y is a vector, only one replicate and variable are 
%               assumed.
%               In the binomial case with local sample sizes M_i,
%               Y is a cell array of length 2, the first of which cantains
%               the matrix or array above containing observed frequencies,
%               and the second of which contains the corresponding local
%               sample sizes.
%  FDPAROBJ ... A functional parameter or fdPar object.  This object 
%               contains the specifications for the functional data
%               object to be estimated by smoothing the data.  See
%               comment lines in function fdPar for details.
%               This argument may also be either a FD object, or a 
%               BASIS object.  If this argument is a basis object, the 
%               smoothing parameter LAMBDA is set to 0.
%
%  Optional arguments are input in pairs:  the first element of the pair
%     is a string specifying the property that the argument value defines,
%     and the second element is the value of the argument
%
%     Valid property/value pairs include:
%
%     Property        Value
%     ----------------------------------------------------------------
%     weight          vector of the same length as the data vector to be
%                     smoothed, containing nonnegative weights to be 
%                     applied to the data values
%     fdnames         A cell array of length 3 with names for
%                       1. argument domain, such as 'Time'
%                       2. replications or cases
%                       3. the function.
%     covariates      A N by Q matrix Z of covariate values used to augment
%                     the smoothing function, where N is the number of
%                     data values to be smoothed and Q is the number of
%                     covariates.  The process of augmenting a smoothing 
%                     function in this way is often called "semi-parametric 
%                     regression".  The default is the empty object [].
%     dfscale         A scalar value multiplying the degrees of freedom
%                     in the definition of the generalized 
%                     cross-validated or GCV criterion for selecting the
%                     bandwidth parameter LAMBDA.  It was recommended by
%                     Chong Gu that this be a number slightly larger than
%                     1.0, such as 1.2, to prevent under-smoothing,
%                     The default is 1.0.
%     family          a character string containing one of:
%                       'normal'  
%                       'binomial'
%                       'poisson'
%                       'gamma'
%                       'inverse gaussian'
%                     the value determines which of the link functions in
%                     the generalized linear model (GLM) family is to be
%                     used.  The default is 'normal'.
%      control        a struct object controlling iterations with members:
%                       epsilon  convergence criterion (default 1e-8)
%                       maxit    max. iterations       (default 25)
%                       trace    output iteration info (0)
%      start          a vector containing starting values for coefficients
%                      
%
%  Returned objects are:
%
%  FDOBJ   ... an object of class fd containing coefficients.
%  DF      ... a degrees of freedom measure.
%  GCV     ... a measure of lack of fit discounted for df.
%                 If the function is univariate, GCV is a vector 
%                 containing the error  sum of squares for each 
%                 function, and if the function is multivariate, 
%                 GCV is a NVAR by NCURVES matrix.
%  SSE     ... the error sums of squares.  
%                 SSE is a vector or matrix of the same size as 
%                 GCV.
%  PENMAT  ... the penalty matrix, if computed, otherwise [].
%  Y2CMAP  ... the matrix mapping the data to the coefficients.
%  ARGVALS ... the input set of argument values.
%  Y       ... the input array containing values of curves

%  Last modified 16 September 2014 by Jim Ramsay

%  check ARGVALS

[argvals, n] = argcheck(argvals);

%  check Y

[y, ncurve, nvar, ndim] = ycheck(y, n);

%  check FDPAROBJ and get FDOBJ and LAMBDA

fdParobj = fdParcheck(fdParobj);
fdobj    = getfd(fdParobj);
lambda   = getlambda(fdParobj);
Lfdobj   = getLfd(fdParobj);

%  check LAMBDA

if lambda < 0, lambda = 0;  end

%  get BASIS and NBASIS

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj) - length(getdropind(basisobj));

%  set default argument values

deffdnames = cell(1,3);
deffdnames{1} = 'arguments';
deffdnames{2} = 'replications';
deffdnames{3} = 'variables';

%  set up the argument names and argument default values

% paramNames = {     'weight' 'fdnames' 'covariates' 'dfscale' 'family'};
% paramDflts = {[]     deffdnames      []        1.0     'normal'};

wtvec      = [];
fdnames    = deffdnames;
covariates = [];
dfscale    = 1;
family     = 'binomial';

% Which style of calling sequence is being used:  
%    name -- value pairs or fixed order?

NameValStyle = true;
if ~isempty(varargin) && nargin <= 7
   va1 = varargin{1};
   if ~ischar(va1) || isempty(va1)
      NameValStyle = false;
   end
end

if NameValStyle
    
    % Process optional number of name -- value pairs.
    
    nargpr = nargin - 3;
    if floor(nargpr/2)*2 ~= nargpr
        error(['The number of argments after the first three ', ...
            'is not an even number.']);
    end
    
    for ipr=5:2:nargin
        ArgName = varargin{ipr-4};
        if     strcmp(ArgName, 'w')         || ...
               strcmp(ArgName, 'wt')        || ...
               strcmp(ArgName, 'wgt')       || ...
               strcmp(ArgName, 'weight')    || ...
               strcmp(ArgName, 'weights')
            wtvec      = varargin{ipr-3};
        elseif strcmp(ArgName, 'f')         || ...
               strcmp(ArgName, 'fdname')    || ...
               strcmp(ArgName, 'fdnames')
            fdnames    = varargin{ipr-3};
        elseif strcmp(ArgName, 'c')         || ...
               strcmp(ArgName, 'cov')       || ...
               strcmp(ArgName, 'covariate') || ...
               strcmp(ArgName, 'covariates')
            covariates = varargin{ipr-3};
        elseif strcmp(ArgName, 'f')         || ...
               strcmp(ArgName, 'fam')       || ...
               strcmp(ArgName, 'family')
            family = varargin{ipr-3};
        elseif strcmp(ArgName, 'd')         || ...
               strcmp(ArgName, 'df')        || ...
               strcmp(ArgName, 'dfscl')     || ...
               strcmp(ArgName, 'dfscale')
            dfscale    = varargin{ipr-3};
        elseif strcmp(ArgName, 'con')       || ...
               strcmp(ArgName, 'control')
            control    = varargin{ipr-3};
        elseif strcmp(ArgName, 's')         || ...
               strcmp(ArgName, 'st')        || ...
               strcmp(ArgName, 'start')
            start    = varargin{ipr-3};
        else
            error('An argument name is unrecognizable.');
        end
    end
else
    
    %  process argument values in fixed order
    
    if nargin >= 4,  wtvec      = varargin{1};   end
    if nargin >= 5,  fdnames    = varargin{2};   end
    if nargin >= 6,  covariates = varargin{3};   end
    if nargin >= 7,  family     = varargin{4};   end
    if nargin == 8,  dfscale    = varargin{5};   end
    if nargin > 7
        error('More than seven non-named arguments found.');
    end
end

%  check WTVEC

[wtvec, onewt] = wtcheck(n, wtvec);
if onewt
    wtvec = ones(n,1);
end

%  check FDNAMES

if ~iscell(fdnames)
    error('smooth_basis_LS:fdnames', ...
          'Optional argument FDNAMES is not a cell array.');
end

if length(fdnames) ~= 3
    error('smooth_basis_LS:fdnames', ...
          'Optional argument FDNAMES is not of length 3.');
end

%  check COVARIATES

q = 0;
if ~isempty(covariates)
    if ~isnumeric(covariates)
        error('smooth_basis_LS:covariates', ...
            'Optional argument COVARIATES is not numeric.');
    end
    if size(covariates,1) ~= n
        error('smooth_basis_LS:covariates', ...
            'Optional argument COVARIATES has incorrect number of rows.');
    end
    q = size(covariates,2);
end

%  ------------------------------------------------------------------
%                set up the linear equations for smoothing
%  ------------------------------------------------------------------

%  set up matrix of basis function values

basismat  = eval_basis(argvals, basisobj);

if n >= nbasis || lambda > 0
    
    %  The following code is for the coefficients completely determined
    
    %  set up additional rows of the least squares problem for the
    %  penalty term.

    basismat0 = basismat;
    y0        = y;
    
    if lambda > 0
%         nderiv  = getnderiv(Lfdobj);
        penmat  = eval_penalty(basisobj, Lfdobj);
        lamRmat = lambda.*penmat;
%         [V,D] = eig(full(penmat));
%         Dvec  = diag(D);
%         [Dsort, Isort] = sort(Dvec, 'descend');
%         Vsort = V(:,Isort);
        %  Check that the lowest eigenvalue in the series that is to be
        %  kept is positive.
%         eiglow = nbasis - nderiv;
%         if Dsort(eiglow) <= 0
%             error('smooth_basis:eig', ...
%                   ['Eigenvalue(NBASIS-NDERIV) of penalty matrix ', ...
%                    'is not positive; check penalty matrix.']);
%         end
%         %  Check that the highest eigenvalue that is not used is small
%         %  relative to the largest eigenvalue.
%         if nderiv > 0 && log10(Dsort(eiglow+1)/Dsort(1)) > -1e-12
%             error('smooth_basis:eig', ...
%                   ['Eigenvalue(NBASIS-NDERIV+1) of penalty matrix ', ...
%                    'is not small relative to eigenvalue(1); ', ...
%                    'check penalty matrix.']);
%         end
        %  Compute the square root of the penalty matrix in the subspace
        %  spanned by the first N - NDERIV eigenvectors
%         ind = 1:eiglow;
%         penfac = Vsort(:,ind)*diag(sqrt(Dsort(ind)));
%         %  Augment basismat by sqrt(lambda).*penfac'
%         basismat = [basismat; sqrt(lambda).*penfac'];
%         %  Augment data vector by n - nderiv 0's
%         if ndim < 3
%             y = [y; zeros(nbasis-nderiv,ncurve)];
%         else
%             for ivar=1:nvar
%                 y(:,:,ivar) = [y(:,:,ivar); zeros(n-nderiv,ncurve)];
%             end
%         end
    else
        lamRmat = [];
    end
    
    %  augment BASISMAT0 and BASISMAT by the covariate matrix 
    %  if it is supplied
    
    if ~isempty(covariates)
        sparsewrd = issparse(basismat0);
        basismat0 = full([basismat0, covariates]);
        basismat  = full([basismat,  covariates]);
        if sparsewrd
            basismat0 = sparse(basismat0);
            basismat  = sparse(basismat);
        end
        if ~isempty(lamRmat)
            lamRmat = [[lamRmat,         zeros(nbasis,q)];
                       [zeros(q,nbasis), zeros(q)       ]];
        end
    end
    
    %  ------------------------------------------------------------------
    %               compute solution using Matlab function glmfit
    %  ------------------------------------------------------------------

    if ndim < 3
        coef  = zeros(nbasis,ncurve);
        dev   = zeros(ncurve,1);
        [coef,dev] = glm_fda(basismat, y, family, lamRmat, wtvec);
    else
        coef  = zeros(nbasis,ncurve,nvar);
        dev   = zeros(ncurve,nvar);
        for ivar=1:nvar
            yi = squeeze(y(:,:,ivar));
            [coefi,devi] = glm_fda(basismat, yi, family, lamRmat, wtvec);
            coef(:,:,ivar) = coefi;
            dev(:,ivar)    = devi;
            stats{:,ivar}  = statsi;
        end
    end
            
    %  compute basismat*R^{-1}
    
    if isempty(lamRmat)
        M = basismat'*basismat;
    else
        M = basismat'*basismat + lamRmat;
    end
    
    %  compute map from y to c
    
    y2cMap = M\basismat';
    
    %  compute degrees of freedom of smooth
    
    df = sum(diag(basismat*y2cMap));
 
else
    error(['The number of basis functions exceeds the number of ', ...
           'points to be smoothed.']);    
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute error sum of squares

if ndim < 3
    yhat = basismat0 * coef;
    SSE  = sum((y0 - yhat).^2);
else
    SSE = zeros(nvar,ncurve);
    for ivar = 1:nvar
        coefi = squeeze(coef(:,:,ivar));
        yhati = basismat * coefi;
        yi    = squeeze(y(:,:,ivar));
        SSE(ivar,:) = sum((yi - yhati).^2);
    end
end

%  compute  GCV index

if df < n
    gcv = (SSE./n)./((n - dfscale*df)/n)^2;
else
    gcv = NaN;
end

%  set up the functional data object

if ndim < 3
    fdobj = fd(coef(1:nbasis,:),   basisobj, fdnames);
else
    fdobj = fd(coef(1:nbasis,:,:), basisobj, fdnames);
end

%  set up the regression coefficient matrix beta

if q > 0
    ind = (nbasis+1):(nbasis+q);
    if ndim < 3
        beta = coef(ind,:);
    else
        beta = coef(ind,:,:);
    end
else
    beta = [];
end






