function quadvals = MakeQuadPoints(knots,nquad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functin MakeQuadPoints
%
% Creates quadrature points and values for using Simpson's rule in a
% Bspline basis. 
%
% INPUTS:
%
% knots  - the knots of the basis
%
% nquad  - number of quadrature points to use
%
% OUTPUT:
% 
% quadvals - a 2-dimensional array containing the quadrature points to use
%            in the first dimension and the corresponding quadrature
%            weights in the second.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

breaks  = unique(knots);
nbreaks = length(breaks);

quadpts1 = linspace(breaks(1),breaks(2),nquad)';
quadwts1 = ones(nquad,1);
quadwts1(2:2:nquad-1) = 4;
quadwts1(3:2:nquad-2) = 2;
quadwts1 = ((breaks(2)-breaks(1))/(nquad-1)).*quadwts1/3;
quadvals = [quadpts1, quadwts1];

for i=3:nbreaks
    quadptsi = linspace(breaks(i-1),breaks(i),nquad)';
    quadwtsi = ones(nquad,1);    
    quadwtsi(3:2:nquad-2) = 2;
    quadwtsi(2:2:nquad-1) = 4;
    quadwtsi = ((breaks(i)-breaks(i-1))/(nquad-1)).*quadwtsi/3;
    quadvals = [quadvals;[quadptsi, quadwtsi]];
end

end