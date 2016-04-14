function coef = getcoef(fdParobj, coef)
%  Retrieve the coefficient matrix from fd object argument.

%  Last modified 4 March 2015

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end
coef = getcoef(fdParobj.fd);