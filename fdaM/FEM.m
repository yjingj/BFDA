function basismat = FEM(evalarg, p, t, nderiv)
%  FEM computes the basis matrix for the FEM basis.  Each row of this 
%  matrix corresponds to a value the corresponding row of EVALARG.  
%  The number of columns of the matrix equals the number of rows of NODES.
%  The value in the position in a row corresponding to the first value in 
%  the same row in EVALARG is set to the value in the second position in
%  the corresponding row of EVALARG.  All other values are zero.  
%  The matrix is returned in sparse matrix format.
%
%  Note:  Unlike spline and other basis functions, the values of FEM basis 
%  functions are either 0 or 1, and, consequently, function FEM probably 
%  has few interesting uses in functional data analysis.   It is included 
%  essentially only for completeness.
%
%  Arguments:
%  EVALARG ... An N by 2 matrix.  In any row, the first value
%                must be an integer, and is the index of the node whose
%                value is to be set to the value in the second value.
%                            
%  P       ... The NBASIS by 2 matrix of vertices of triangles containing
%                 the X- and Y-coordinates of the vertices.
%              P may also be provided as a 2 by NBASIS matrix.
%  T       ... The 4 by no. of triangles matrix specifying triangles and 
%              their properties:
%              Columns 1 to 3:   the indices in P of the vertices of each 
%              triangle in counter-clockwise order.
%              Column  4:        the subdomain number of the triangle
%              T may also be provided with 3 columns, in which case the
%              final column is set to ones.  It may also be provided as
%              either a 3 or 4 by number of triangles matrix.
%  NDERIV  ... An order of derivative to be evaluated, with may be
%              0, 1 or 2.

%  Returns:  
%  BASIS   ... An N by NBASIS matrix of FEM basis function values in
%              sparse matrix format.

%  Last modified 14 June 2010 by Jim Ramsay

if nargin < 4,  nderiv = 0;  end

%  check EVALARG

if ~isnumeric(evalarg)
    error('EVALARG is not a numerical object.');
end
    
argdim = size(evalarg);

if length(argdim) ~= 2
    error('EVALARG is not a matrix.');
else
    N = argdim(1);
    m = argdim(2);
end

if N ~= 2 && m ~= 2
    error('EVALARG does not have either two rows or two columns.');
end

if N == 2 && m ~= 2
    %  transpose if EVALARG is two rows and N columns
    evalarg = evalarg';
    N = m;
end

if nderiv ~= 0
    error('NDERIV is not zero.');
end

ind = evalarg(:,1);

if any(ind ~= floor(ind))
    error('Non-integer values in the first column of EVALARG');
end

%  check dimensions of P, E and T and transpose if necessary

if     size(t,2) == 3
    t = [t, ones(size(t,1),1)];
elseif size(t,1) == 3
    t = [t; ones(size(t,1),1)'];
else
    error('Neither dimension of T is either 3 or 4.');
end

if size(p,1) ~= 2 && size(p,2) ~= 2 || ...
   size(t,1) ~= 4 && size(t,2) ~= 4 
   error('Dimensions of at least one of P and T are not correct.');
end

if size(p,2) ~= 2 && size(p,1) == 2
    p = p';
end
if size(t,2) ~= 4 && size(t,1) == 4
    t = t';
end

%  set up the nodes and node mesh

nodes    = makenodes(p,t);

nbasis   = size(nodes,1);

basismat = sparse(N,nbasis);

%  check indices

errorwrd = 0;
for i=1:N
    if ind(i) < 1 || ind(i) > nbasis
        errorwrd = 1;
    end
end
if errorwrd
    error('One or more indices in EVALARG are out of range.');
end

for i=1:N
    basismat(i,ind(i)) = evalarg(i,2);
end


