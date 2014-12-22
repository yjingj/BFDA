%Generate 2 eigenfunctions
function [eigen] = genEigenFun_1(t,der)
     eigen = zeros(2,length(t));
     if nargin == 1|| der == 0
       eigen(1,:) = sqrt(2)*sin(pi*t);
       eigen(2,:) = sqrt(2)*cos(pi*t);
     elseif nargin == 2

        if der == 1
           eigen(1,:) = pi*sqrt(2)*cos(pi*t);
           eigen(2,:) = -sqrt(2)*pi*sin(pi*t);
        elseif der == 2
           eigen(1,:) = -sqrt(2)*pi^2*sin(pi*t);
           eigen(2,:) = -sqrt(2)*pi^2*cos(pi*t);
        end           

     end    

end
