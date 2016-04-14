function bifdParobj = bifdParcheck(bifdParobj)

if ~isa_bifdPar(bifdParobj) 
    if isa_bifd(bifdParobj) || isa_basis(bifdParobj)
        bifdParobj = bifdPar(bifdParobj);
    else
        error(  ...
        ['BIFDPAROBJ is not a bivariate functional parameter object, ', ...
         'not a bivariate functional data object, and ', ...
         'not a bivariate basis object.']);
    end
end
