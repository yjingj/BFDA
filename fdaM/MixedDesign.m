function Z = MixedDesign(a, b, r)
%  MixedDesign sets up the design matrix for one between and one within
%  factor in an ANOVA.
%  Arguments:
%  a  ... number of levels of the within  factor A
%  b  ... number of levels of the between factor B
%  r  ... vector of length b containing the numbers of subjects
%         in each level of factor B

%  Last modified 6 Feb. 2007

N = a*sum(r);            %  number of dependent variable values
p = 1 + a + b + a*b + N; %  number of parameters
M = 3 + b;               %  number of zero sum constraints