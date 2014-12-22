%function [newz] = mapX2d(x,y,z,newx,newy)
%Map (x,y,z) to (newx,newy,newz) 
%i) x,y,z can be all 1*n vectors
%ii)x : 1 * p vector
%   y : 1 * q vector
%   z(i,j) corresponds to x(i) and y(j)
%   where i = 1,...,p; j = 1,...,q
%   newx: 1 * m 
%   newy: 1 * m
%   newz: 1 * m
%  example
%  i)
%   x = 1:3;
%   y = 2:4;
%   z = 0.1:0.1:0.3;
%   mapX2d(x,y,z,[1 3],[2 4])
%
%  ii)
%  x = 1:3;
%  y = 2:4;
%  z = [0.1:0.1:0.3; 0.2:0.1:0.4;0.3:0.1:0.5];
%  newx = [1 1 2 3];
%  newy = [2 3 4 2];
%  mapX2d(x,y,z,newx,newy)
function [newz] = mapX2d(x,y,z, newx,newy)

     if size(z,1) > 1 && size(z,2) > 1
         [xx,yy] = meshgrid(x,y);             
         xx = xx(:)';
         yy = yy(:)';
         zz = z';
         zz = zz(:)';
         [newz] = mapX2d(xx,yy,zz,newx,newy);

     else
       
       newz = zeros(1,length(newx));
       for i = 1:length(newx)
         id = find(x == newx(i) & y == newy(i));
         if isempty(id)
             newz(i) = NaN;
         else
	     newz(i) = z(id);
         end
       end

     end

end
