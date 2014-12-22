function [res]=romb2(x,y, decdigs,varargin)
%  [res]=romb2(x,y, decdigs)
%  Use 2^(decidigs) trapezoids to estimate the integral
%  default value is 10, where the original function 
%  works for a user-defined function from a to b
%  Input 
%  x    :    1 x n vector of time grids
%  y    :    1 x n vector of function value that evaluated at x
%            or n x m matrix, each of y(:,j), j = 1,...,m is evaluated at x 
%  decdigs:  positive integer
%            2^decdigs number of trapezoids will be used in
%            the integration, default is 10    
%  Output 
%  res:      real number of the numerical integral corresponding
%            to the input.  
%            or m x 1 vector
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
  m = size(y,2);
  res = zeros(m,1);
  for i = 1:m
    res(i) = romb(x,y(:,i),decdigs);
  end
end
 
