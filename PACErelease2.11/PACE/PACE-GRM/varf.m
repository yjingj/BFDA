% variance function
function v = varf(mux,family,param)

    if strcmp(family,'binomial') == 1
        mux(mux<=0) = exp(-5);
        mux(mux>=1) = 1-exp(-5);
        v = mux.*(1-mux);
    elseif strcmp(family,'poisson') == 1
        v = mux;
    end
