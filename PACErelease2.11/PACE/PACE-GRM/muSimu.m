function mut = muSimu(t,T)

    id = find(t < 0 | t > T);
    if ~isempty(id)
       fprintf(1,['Warning: t must be in [0, ' num2str(T) ']. Invalid t''s are removed!\n']);
       t(id) = [];
    end

    mut =  2*sqrt(2/T)*sin(2*pi*t/T);
