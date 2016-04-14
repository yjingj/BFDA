function bifdobj = getbifd(bifdParobj)
%GETBIFD   Extracts the bivariate functional data object from 
%   a bivariate functional parameter object.

%  last modified 23 November 2012

if ~isa_bifdPar(bifdParobj)
    error('Argument is not a bivariate functional parameter object.');
end

bifdobj = bifdParobj.bifd;


