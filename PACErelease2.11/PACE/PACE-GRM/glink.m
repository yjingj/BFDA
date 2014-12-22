% link function
function gfun = glink(x,family,param)

    if strcmp(family,'binomial') == 1
        gfun = exp(x)./(1+exp(x));
    elseif strcmp(family,'poisson') == 1
        gfun = exp(x);
    end


