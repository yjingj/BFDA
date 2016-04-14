function K0 = stiff0(nodeStruct)

%STIFF0 produces the stiffness matrix for the equation which
%ensures that del(u)=w in the mixed method solution to the 
%biharmonic (or minimal curvature) problem, solved by quadratic triangles.
%
%Input: FMESH is the triangularization matrix produced by
%	calling fmesh=timtri(points).  Note that points is
%	an nver*2 matrix of points (x,y).
%       PSI is the matrix of shape function derivatives
%       produced by the call qder.
%
%Output: K0: the output is an nnod*nnod matrix out which is
%        the sum of the nele element 0-stiffness matrices.
%        For each element i and each 1<=j,k<=6, the integral of the product 
%        of the j'th and k'th shape functions over the i'th element is
%        computed.  Then that value is the (fmesh(i,j),fmesh(i,k))'th entry
%        of the i'th elemental stiffness matrix.
% Last modified on September 11, 2000.

%   make matrix psi representing basis functions as columns.
%   The 6-tuple (a,b,c,d,e,f) represents the function
%	 a + b*L2 + c*L3 + d*(L2)^2 + e*L2*L3 + f*(L3)^2, where L2=x,
%	 L3=y, L1=1-L2-L3 are the barycentric coordinates of the universal
%	 triangular element.

psi1 = [ 1 ;-3 ;-3 ; 2 ; 4 ; 4 ];
psi2 = [ 0 ; 4 ; 0 ;-4 ;-4 ; 0 ];
psi3 = [ 0 ;-1 ; 0 ; 2 ; 0 ; 0 ];
psi4 = [ 0 ; 0 ; 0 ; 0 ; 4 ; 0 ];
psi5 = [ 0 ; 0 ;-1 ; 0 ; 0 ; 2 ];
psi6 = [ 0 ; 0 ; 4 ; 0 ;-4 ;-4 ];
psi  = [ psi1 psi2 psi3 psi4 psi5 psi6 ];

order     = nodeStruct.order;
nodes     = nodeStruct.nodes;
nodeindex = nodeStruct.nodeindex;
Jvec      = nodeStruct.J;

nele   = size(nodeindex,1);
nnod   = size(nodes,1);
K0     = sparse(nnod,nnod);

w      = diag([9/80 repmat((155+sqrt(15))/2400,1,3) ...
                 repmat((155-sqrt(15))/2400,1,3)]);
x      = 1/21*[7 6+sqrt(15) 9-2*sqrt(15) 6+sqrt(15) ...
                  6-sqrt(15) 9+2*sqrt(15) 6-sqrt(15)]';
y      = 1/21*[7 6+sqrt(15) 6+sqrt(15) 9-2*sqrt(15)  ...
                  6-sqrt(15) 6-sqrt(15) 9+2*sqrt(15)]';
q      = [repmat(1,7,1) x y x.*x x.*y y.*y];
Khat   = psi'*q'*w*q*psi;

J      = (nodes(nodeindex(:,3),1)-nodes(nodeindex(:,1),1)).*...
         (nodes(nodeindex(:,5),2)-nodes(nodeindex(:,1),2))-...
	     (nodes(nodeindex(:,5),1)-nodes(nodeindex(:,1),1)).*...
	     (nodes(nodeindex(:,3),2)-nodes(nodeindex(:,1),2));

for el=1:nele
    ind = nodeindex(el,:);
    K0(ind,ind) = K0(ind,ind) +	Khat*J(el);
end;


