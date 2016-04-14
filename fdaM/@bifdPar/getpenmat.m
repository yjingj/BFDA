function [penmats, penmatt] = getpenmat(bifdParobj)
%GETPENMAT   Extracts the penmat value from 
%   a functional parameter object.

%  last modified 26 November 2012

if ~isa_bifdPar(bifdParobj)
    error('Argument is not a bivariate functional parameter object');
end
penmats = bifdParobj.penmats;
penmatt = bifdParobj.penmatt;


