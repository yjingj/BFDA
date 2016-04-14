function bwtcell = getbwtcell(Lfdobj)
%  GETWFDCELL   Extracts the weight function cell object from LFDOBJ.

%  last modified 18 June 2013

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

bwtcell = Lfdobj.bwtcell;


