function [res]=romb(x,y, decdigs,varargin)
%  [res]=romb(x,y, decdigs)
%  Use 2^(decidigs) trapezoids to estimate the integral
%  default value is 10, where the original function 
%  works for a user-defined function from a to b
%  Input 
%  x    :    1 x n vector of time grids
%  y    :    1 x n vector of function value that evaluated at x
%  decdigs:  positive integer
%            2^decdigs number of trapezoids will be used in
%            the integration, default is 10    
%  Output 
%  res:      real number of the numerical integral corresponding
%            to the input.  
%  ex: 
%   x = 0:0.01:pi;
%   y = sin(x);
%   res = romb(x,y);
%   or 
%   res = romb(x,y,8);
%
%  Modified by Bitao from the following source:
% 
%ROMBINT	 Numerical evaluation of an integral using the Romberg method.
%
%   ---------------------------------------------------------
%   Author: Martin Kacenak,
%           Department of Applied Informatics and Process Control,
%           Faculty of BERG, Technical University of Kosice,
%           B.Nemcovej 3, 04200 Kosice, Slovak Republic
%   E-mail: ma.kac@post.cz
%   Date:   posted February 2001, updated June 2006.
%
if nargin<3, decdigs=10; end
if nargin<2
   warning ('Error in input format')
else
   rom=zeros(2,decdigs);
   romall=zeros(1,(2^(decdigs-1))+1);   
%   romall=feval(funfcn,a:(b-a)/2^(decdigs-1):b,varargin{:});
   a = min(x);
   b = max(x);
   romall = interp1(x,y, a:(b-a)/2^(decdigs-1):b, 'spline');
   h=b-a;
   rom(1,1)=h*(romall(1)+romall(end))/2;
   for i=2:decdigs
      st=2^(decdigs-i+1);
      rom(2,1)=(rom(1,1)+h*sum(romall((st/2)+1:st:2^(decdigs-1))))/2;
      for k=1:i-1
         rom(2,k+1)=((4^k)*rom(2,k)-rom(1,k))/((4^k)-1);
      end
      rom(1,1:i)=rom(2,1:i);
      h=h/2;
   end
   res=rom(1,decdigs);
end
 
