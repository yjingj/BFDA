function [b,R] = wfit(y,x,sw)
% Perform a weighted least squares fit
p  = size(x,2);
yw = y .* sw;
xw = x .* sw(:,ones(1,p));
% No pivoting, no basic solution.  We've removed dependent cols from x, and
% checked the weights, so xw should be full rank.
[Q,R] = qr(xw,0);
b = R \ (Q'*yw);
