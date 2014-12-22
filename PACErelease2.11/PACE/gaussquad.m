%[res] = gaussquad(x,y,n)
%Compute the integral through Gaussian Quadrature.
%It first calls the lgwt.m written by Greg von Winckel
%to obtain the Legendre-Gauss weights and nodes for computing
%the definite integral of a continuous function on some interval
%[a,b], here, a = min(x) and b = max(x). The realization values
%for the function is y.
%Input:
%x:  1 * n vector 
%y:  n * k vector, contains function values that evaluated at x
%n:  integer denotes the nth Legendre polynomials 
%Output:
% res : k * 1 vector of resulting integrated value    
%ex: x = 0:0.01:pi;
%    y = [sin(x)' cos(x)' sin(x)'];
%    res = gaussquad(x,y)
%    or 
%    res = gaussquad(x,y,5)
function [res] = gaussquad(x,y,n)

     if nargin < 3
        n = 30;
     end
     a = min(x);
     b = max(x);
     [x_node,w] = lgwt(n,a,b);
     y_node = interp1(x,y, x_node, 'spline');
     res = y_node'*w;

end

