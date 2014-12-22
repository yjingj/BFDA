%function [d] = getDist(x1,x2,y1,y2, toSort)
%This function obtain the Euclidean distance between (x1,y1) and (x2,y2),
%it can be either a point or vector of points.
%Input arguments:
%x1,y1,x2,y2: 1 * n vectors
%toSort: 1   sort the resulting distance
%         0   Do not sort the resulting distance
%Output arguments: 
%d : 1 * n vector of the resulting Euclidean distances between (x1,y1) and (x2,y2)
function [d] = getDist(x1,x2, y1, y2, toSort)

     d = sqrt((x1-x2).^2+(y1-y2).^2);
     if toSort == 1
       d = sort(d);
     end

end
