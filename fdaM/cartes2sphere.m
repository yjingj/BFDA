function [angle, radius] = cartes2sphere(coord)
%  CARTES2SPHERE converts Cartesian coordinates of a point on a sphere
%  in a space of dimension n + 1 into n angles and the radius

%  Last modified 11 August 2012 by Jim Ramsay

n = length(coord) - 1;
np1 = n + 1;
coordsqr = coord.^2;
ind = np1:-1:1;
rtcsum = sqrt(cumsum(coordsqr(ind)));
angle = zeros(n,1);
angle(1:n-1) = acot(coord(1:n-1)./rtcsum(ind(2:n)));
zeroind = find(coord(1:n-1) == 0 & rtcsum(ind(2:n)) == 0);
if coord(min(zeroind)-1) > 0
    angle(zeroind) = 0;
else
    angle(zeroind) = pi;
end    
if angle(n-1) == 0
    angle(n) = 0;
elseif angle(n-1) == pi
    angle(n) = pi;
else
    angle(n) = 2*acot((rtcsum(2) + coord(n))/coord(n+1));
end
radius = rtcsum(n+1);
negind = angle < 0;
angle(negind) = angle(negind) + pi;
