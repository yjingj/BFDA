%function [newy] = mapX1d(x,y,newx)
%Map (x,y) to (newx,newy)
%x, y : vectors of 1 * n
%newx : vector of 1 * m
%newy : vector of 1*m
function [newy] = mapX1d(x,y,newx)

     [ignore, loc] = ismember(newx, x);
     if size(y,1) == 1
        newy = y(loc);
     elseif size(y,1) > 1
        newy = y(loc,:);
     else
        fprintf(1,'y cannot be empty!\n');
     end   
end
