function fdParobjnew = putcoef(fdParobj, coef)
%  Replace the fd object coefficent matrix by argument COEF.

%  Last modified 4 March 2015

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end
fdParobjnew = fdParobj;
basisobj = getbasis(fdParobj.fd);
nbasis   = size(getcoef(fdParobj.fd),1);
if size(coef,1) == nbasis
    fdParobjnew.fd = fd(coef, basisobj);
else
    error(['New coefficient matrix has number of rows not equal to ', ...
           'that of old coefficient matrix.']);
end
