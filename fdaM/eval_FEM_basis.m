function evalmat = eval_FEM_basis(Xvec, Yvec, basisobj, nderivs)
% EVAL_FEM_BASIS evaluates the FEM basis object at points in (Xvec,Yvec)
%
% Arguments:
% Xvec     ... an vector of x-coordinates of length N.  This may also be
%              a matrix with two columns, in which case the first column
%              is assigned to Xvec and the second to Yvec.
% Yvec     ... an vector of y-coordinates of length N.
% BASISOBJ ... A functional basis object of the FEM type.
% NDERIVS  ... A vector of length 2 containing the orders of derivatives
%              for X and Y, respectively.
% 
%        output:
% EVALMAT   an array of the same size as Xvec and Yvec containing the value  
%           of FELSPLOBJ at (Xvec,Yvec).
%
%  Last modified on 29 October 2014 by Jim ramsay.

%  Set up the arguments if the first argument is a matrix with two
%  columns

%  set default values

if nargin < 4,  nderivs = zeros(1,2);  end

%  allow for matrix input of point coordinates

if size(Xvec,2) == 2
    if nargin < 2
        error(['First argument is a coordinate matrix and ', ...
               'the second argument is not supplied.']);
    end
    basisobj = Yvec;
    Yvec     = Xvec(:,2);
    Xvec     = Xvec(:,1);
else
    if nargin < 3
        error(['First and second arguments are coordinate vectors ', ...
               'and the third argument is not supplied.']);
    end
end

%  check the type of BASISOBJ

if ~strcmp(getbasistype(basisobj), 'FEM')
    error('The basis object for BASISOBJ is not of type FEM.');
end
nbasis = getnbasis(basisobj);

%  check Xvec

if ~isa(Xvec,'double')
   error('Xvec is not a numerical array')
else
   Xvec  = Xvec(:);     % treat Xvec as a column vector
end

%  check Yvec

if ~isa(Yvec,'double')
   error('Yvec is not a numerical array')
elseif length(Yvec(:))~=length(Xvec)
   error('Yvec is not the same length as Xvec')
else
   Yvec=Yvec(:);     % treat Yvec as a column vector
end

%  check derivatives

if length(nderivs) ~= 2
    error('NDERIVS not of length 2.');
end
if sum(nderivs) > 2
    error('Maximum derivative order is greater than two.');
end

N = length(Xvec);

%  Augment Xvec and Yvec by one's for computing barycentric coordinates 

Pgpts = [ones(N,1) Xvec Yvec];

%  get nodes and index

params    = getbasispar(basisobj);
p         = params.p;
t         = params.t;
t         = t(:,1:3);
order     = params.order;
nodes     = params.nodes;
nodeindex = params.nodeindex;
Jvec      = params.J;

% 1st, 2nd, and 3rd vertices of triangles

if     order == 2
    v1 = nodes(nodeindex(:,1),:);
    v2 = nodes(nodeindex(:,3),:);
    v3 = nodes(nodeindex(:,5),:);
elseif order == 1
    v1 = nodes(nodeindex(:,1),:);
    v2 = nodes(nodeindex(:,2),:);
    v3 = nodes(nodeindex(:,3),:);
else
    error('ORDER is neither 1 nor 2.');
end

% denominator of change-of-coordinates change matrix

modJac    = Jvec;
ones3     = ones(3,1);
modJacMat = modJac*ones3';

% 1st, 2nd, and 3rd columns of transformations to barycentric coordinates,
% with a row for each vertex

M1 = [ v2(:,1).*v3(:,2) - v3(:,1).*v2(:,2) ...
       v2(:,2)-v3(:,2)                     ...
       v3(:,1)-v2(:,1)]./modJacMat/2;
M2 = [ v3(:,1).*v1(:,2) - v1(:,1).*v3(:,2) ...
       v3(:,2)-v1(:,2)                     ...
       v1(:,1)-v3(:,1)]./modJacMat/2;
M3 = [ v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2) ...
       v1(:,2)-v2(:,2)                     ...
       v2(:,1)-v1(:,1)]./modJacMat/2;

% Identify triangles containing points in vector (Xvec(i),Yvec(i))
% if no triangle contains a point, ind(i) is NaN

tricoef = tricoefCal(p, t);
ind = zeros(N,1);
for i=1:N
    ind(i) = insideIndex(Xvec(i), Yvec(i), p, t, tricoef);
end

%  interpolate values

evalmat = sparse([], [], [], N, nbasis, 0);

for i=1:N
    indi = ind(i);
    if ~isnan(indi)
        %  change to barycentric coordinates
        baryc1 = (M1(indi,:).*Pgpts(i,:))*ones3;
        baryc2 = (M2(indi,:).*Pgpts(i,:))*ones3;
        baryc3 = (M3(indi,:).*Pgpts(i,:))*ones3;
        if order == 2
            if    sum(nderivs) == 0
                evalmat(i,nodeindex(indi,1)) = 2*baryc1.^2 - baryc1;
                evalmat(i,nodeindex(indi,2)) = 2*baryc2.^2 - baryc2;
                evalmat(i,nodeindex(indi,3)) = 2*baryc3.^2 - baryc3;
                evalmat(i,nodeindex(indi,4)) = 4*baryc1.* baryc2;
                evalmat(i,nodeindex(indi,5)) = 4*baryc2.* baryc3;
                evalmat(i,nodeindex(indi,6)) = 4*baryc3.* baryc1;
            elseif nderivs(1) == 1 && nderivs(2) == 0
                evalmat(i,nodeindex(indi,1)) = (4*baryc1 - 1).*M1(indi,2);
                evalmat(i,nodeindex(indi,2)) = (4*baryc2 - 1).*M2(indi,2);
                evalmat(i,nodeindex(indi,3)) = (4*baryc3 - 1).*M3(indi,2);
                evalmat(i,nodeindex(indi,4)) = 4*baryc2.*M1(indi,2) + ...
                                               4*baryc1.*M2(indi,2);
                evalmat(i,nodeindex(indi,5)) = 4*baryc3.*M2(indi,2) + ...
                                               4*baryc2.*M3(indi,2);
                evalmat(i,nodeindex(indi,6)) = 4*baryc1.*M3(indi,2) + ...
                                               4*baryc3.*M1(indi,2);
            elseif nderivs(1) == 0 && nderivs(2) == 1
                evalmat(i,nodeindex(indi,1)) = (4*baryc1 - 1).*M1(indi,3);
                evalmat(i,nodeindex(indi,2)) = (4*baryc2 - 1).*M2(indi,3);
                evalmat(i,nodeindex(indi,3)) = (4*baryc3 - 1).*M3(indi,3);
                evalmat(i,nodeindex(indi,4)) = 4*baryc2.*M1(indi,3) + ...
                                               4*baryc1.*M2(indi,3);
                evalmat(i,nodeindex(indi,5)) = 4*baryc3.*M2(indi,3) + ...
                                               4*baryc2.*M3(indi,3);
                evalmat(i,nodeindex(indi,6)) = 4*baryc1.*M3(indi,3) + ...
                                               4*baryc3.*M1(indi,3);
            elseif nderivs(1) == 1 && nderivs(2) == 1
                evalmat(i,nodeindex(indi,1)) = 4*M1(indi,2)*M1(indi,3);
                evalmat(i,nodeindex(indi,2)) = 4*M2(indi,2)*M2(indi,3);
                evalmat(i,nodeindex(indi,3)) = 4*M3(indi,2)*M3(indi,3);
                evalmat(i,nodeindex(indi,4)) = 4*M2(indi,2)*M1(indi,3)+ ...
                                               4*M2(indi,3)*M1(indi,2);
                evalmat(i,nodeindex(indi,5)) = 4*M3(indi,2)*M2(indi,3)+ ...
                                               4*M3(indi,3)*M2(indi,2);
                evalmat(i,nodeindex(indi,6)) = 4*M1(indi,2)*M3(indi,3)+ ...
                                               4*M1(indi,3)*M3(indi,2);
            elseif nderivs(1) == 2 && nderivs(2) == 0
                evalmat(i,nodeindex(indi,1)) = 4*M1(indi,2)*M1(indi,2);
                evalmat(i,nodeindex(indi,2)) = 4*M2(indi,2)*M2(indi,2);
                evalmat(i,nodeindex(indi,3)) = 4*M3(indi,2)*M3(indi,2);
                evalmat(i,nodeindex(indi,4)) = 8*M2(indi,2)*M1(indi,2);
                evalmat(i,nodeindex(indi,5)) = 8*M3(indi,2)*M2(indi,2);
                evalmat(i,nodeindex(indi,6)) = 8*M1(indi,2)*M3(indi,2);
            elseif nderivs(1) == 0 && nderivs(2) == 2
                evalmat(i,nodeindex(indi,1)) = 4*M1(indi,3)*M1(indi,3);
                evalmat(i,nodeindex(indi,2)) = 4*M2(indi,3)*M2(indi,3);
                evalmat(i,nodeindex(indi,3)) = 4*M3(indi,3)*M3(indi,3);
                evalmat(i,nodeindex(indi,4)) = 8*M2(indi,3)*M1(indi,3);
                evalmat(i,nodeindex(indi,5)) = 8*M3(indi,3)*M2(indi,3);
                evalmat(i,nodeindex(indi,6)) = 8*M1(indi,3)*M3(indi,3);
            else
                error('Inadmissible derivative orders.');
            end
        else
            if sum(nderivs) == 0
                evalmat(i,nodeindex(indi,1)) = baryc1;
                evalmat(i,nodeindex(indi,2)) = baryc2;
                evalmat(i,nodeindex(indi,3)) = baryc3;
            elseif nderivs(1) == 1 && nderivs(2) == 0
                evalmat(i,nodeindex(indi,1)) = M1(indi(1),2);
                evalmat(i,nodeindex(indi,2)) = M1(indi(2),2);
                evalmat(i,nodeindex(indi,3)) = M1(indi(3),2);
            elseif nderivs(1) == 0 && nderivs(2) == 1
                evalmat(i,nodeindex(indi,1)) = M1(indi(1),3);
                evalmat(i,nodeindex(indi,2)) = M1(indi(2),3);
                evalmat(i,nodeindex(indi,3)) = M1(indi(3),3);
            else
                error('Inadmissible derivative orders.');
            end
        end
    end
end

