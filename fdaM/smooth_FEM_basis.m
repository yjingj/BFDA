function [smooth_fd, df, gcv, beta, SSE, laplace_fd, K0mat, K1mat, ...
          y2cMap, Xvec, Yvec, data] = ...
               smooth_FEM_basis(Xvec, Yvec, data, basisobj, lambda, ...
                                wtvec, covariates)
% SMOOTH_FEM_FD Smooths spatial data observed at a set of points, 
%  possibly using covariate values as additional contributors to the fit to 
%  the data values.  The basis object is of type FEM defined by a 
%  triangular mesh within one or more bounding polygons.  The basis 
%  functions are either piecewise linear or piecewise quadratic functions.
%  The piecewise linear functions are nonzero only over the hexagon 
%  surrounding a VERTEX and with value one at the VERTEX.  The piecewise
%  quadratic functions include these three functions, and in addition
%  including the additional three functions that reach one at the midpoints
%  of the edges of triangles.  The points at which basis functions reach 
%  are called NODES of the mesh. Thus, vertices and nodes correspond for
%  piecewise linear basis functions, but the nodes also include midpoints
%  in the case of piecewise quadratidc basis functions.  Note that all
%  basis functions are continuous across edges and vertices of triangles,
%  but only differentiable within a triangle.  In the piecewise quadratic
%  case, basis functions also so have nonzero second partial derivatives
%  within triangles, but not along edges or at vertices.
%
%  The roughness penalty on the fitting function u(x,y) uses the Laplacian 
%  operator 
%             di^2 u/di x^2 + di^2 u/ di y^2
%  which measures the local curvature at point (x,y).

%  The coordinates of the observation points are entered in vectors XVEC
%  and YVEC, and may be at locations that do not correspond to any node
%  in the mesh.
%
%  Arguments:
% XVEC       ... an vector of x-coordinates of length N.  This may also be
%                a matrix with two columns, in which case the first column
%                is assigned to Xvec and the second to Yvec.
% YVEC       ... an vector of y-coordinates of length N.
% DATA       ... An N by NSURF data matrix containing observations at N
%                points for each of NSURF surfaces.
% BASISOBJ   ... a functional basis object having a basis of the FEM type.
% LAMBDA     ... a scalar smoothing parameter, MUST be positive
% WTVEC      ... vector of length N, containing nonnegative weights to be  
%                applied to the data values, or a symmetric positive  
%                definite matrix of order N.
% COVARIATES ... a design matrix with rows corresponding to observations
%                points and columns corresponding to covariates.
%
%     Output:
% SMOOTH_FD  ...  A FD object of the FEM type defined by the coefficient
%                 vector resulting from smoothing
% LAPLACE_FD ...  A FD object of the FEM type for the value of the 
%                 Laplace operator
%

% Last modified on 1 November 2014 by Jim Ramsay

%  assign defaults to missing arguments

if nargin < 7,  covariates = [];     end
if nargin < 6,  wtvec      = [];     end
if nargin < 5,  lambda     = 1e-12;  end

if nargin < 4
    error('Less than four arguments supplied.');
end

%  ---------------------------------------------------------------
%                      Check arguments
%  ---------------------------------------------------------------

[n,nsurf] = size(data);

%  check XVEC, YVEC and DATA

if length(Xvec) ~= n || size(Xvec,2) ~= 1
    error('XVEC is not a column vector of n X-coordinate values.');
end
if length(Yvec) ~= n || size(Yvec,2) ~= 1
    error('YVEC is not a column vector of n X-coordinate values.');
end

%  Check BASISOBJ

if ~isa_basis(basisobj)
   error('BASISOBJ is not a functional basis object.');
end
if ~strcmp(getbasistype(basisobj), 'FEM')
    error('The basis object for BASISOBJ is not of type FEM.');
end

%  LAMBDA

if  ~isa(lambda,'double')
   error('LAMBDA is not numeric')
elseif size(lambda) ~= [1 2]
   error('LAMBDA is not a scalar')
end

%  check WTVEC

[wtvec, onewt, matwt] = wtcheck(n, wtvec);

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

%  ---------------------------------------------------------------
% construct projection matrix on the space spanned by the columns of the 
% design matrix DESMAT and remove projection of data on DESMAT
% from DATA
%  ---------------------------------------------------------------

if matwt
    wtmat  = wtvec;
    wtfac  = chol(wtmat);
else
    wtmat  = diag(wtvec);
    wtfac  = sqrt(wtmat);
end

if ~isempty(covariates)
    [Q,R] = qr(wtfac*covariates, 0);
    beta  = R\(Q'*(wtfac*data));
    H     = Q * Q';
    data  = data - Q*Q'*data;
else
    beta = [];
end

%  -------------------------------------------------------------------
%     Set up n by nbasis matrix of basis values at observation points
%  -------------------------------------------------------------------

basismat = eval_FEM_basis(Xvec, Yvec, basisobj);

%  ------------------------------------------------------------------
%                Set up the linear equations for smoothing
%  ------------------------------------------------------------------

params   = getbasispar(basisobj);
numnodes = size(params.nodes,1);
indnodes = 1:numnodes;

%  extract quantities required for setting up mass and stiffness matrices

nodeStruct.order     = params.order;
nodeStruct.nodes     = params.nodes;
nodeStruct.nodeindex = params.nodeindex;
nodeStruct.J         = params.J;
nodeStruct.metric    = params.metric;

%  ---------------------------------------------------------------
%                 construct mass matrix K0 
%  ---------------------------------------------------------------

K0mat = mass(nodeStruct);

%  ---------------------------------------------------------------
%              construct stiffness matrix K1
%  ---------------------------------------------------------------

K1mat = stiff1(nodeStruct);

%  ---------------------------------------------------------------
% construct the block diagonal matrix L, having upper left block 
% given by I-H and zero otherwise
%  ---------------------------------------------------------------

%  ---------------------------------------------------------------
% construct vector b for system Ax = b
%  ---------------------------------------------------------------
    
Bmat = sparse(zeros(numnodes*2,nsurf));
Bmat(indnodes,:) = basismat'*wtmat*data;

%  ---------------------------------------------------------------
% construct matrix A for system Ax=b.
%  ---------------------------------------------------------------

Lmat = basismat'*wtmat*basismat;
Amat = [ Lmat    -lambda*K1mat; ...
         K1mat           K0mat ];

% solve system

coefmat = Amat\Bmat;  

%  ---------------------------------------------------------------
%   Make output functional data objects SMOOTH_FD and LAPLACE_FD
%  ---------------------------------------------------------------
    
coef1 = coefmat(indnodes,:);
coef2 = coefmat(indnodes+numnodes,:);
    
smooth_fd  = fd(coef1, basisobj);

if nargout > 1
    laplace_fd = fd(coef2, basisobj);
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute map from y to c

Mmat = sparse(zeros(numnodes*2,n));
Mmat(indnodes,:) = basismat'*wtmat;
temp   = Amat\Mmat;
y2cMap = temp(indnodes,:);

%  compute degrees of freedom of smooth

df = trace(y2cMap*basismat);

%  compute error sum of squares

datahat = basismat * coef1;

SSE = sum(wtfac*(data - datahat).^2);

%  compute  GCV index

if df < n
    gcv = (SSE./n)./((n - df)/n)^2;
else
    gcv = NaN;
end

