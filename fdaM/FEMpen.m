function penaltymat = FEMpen(basisobj, Lfdobj)

if nargin < 2, Lfdobj = int2Lfd(1);  end

if ~isinteger(Lfdobj)
    error('Argument Lfdobj must be an integer.');
end

nderiv = getnderiv(Lfdobj);

if nderiv < 0 || nderiv > 2
    error('Value of argument Lfdobj is not 0, 1 or 2.');
end

if ~strcmp(class(basisobj), 'basis')
    error('Argument basisobj is not a basis object.');
end

if ~strcmp(getbasistype(basisobj), 'FEM')
    error('Argument basisobj is not of type FEM.');
end

params = getbasispar(basisobj);
p = params.p;
t = params.t;

nbasis = getnbasis(basisobj);
np = size(p,2);

%  check order of basis:  must be either 1 or 2

if nbasis == np
    order = 1;
    if nderiv == 2
        error('Second derivative cannot be computed for order 1 basis.');
    end
elseif nbasis == 2*np
    order = 2;
else
    error('NBASIS inconsistent with order being 1 or 2.');
end

if nderiv > 0
    penaltymat = stiff(p, t, order, nderiv);
else
    penaltymat = mass(p, t, order);
end
