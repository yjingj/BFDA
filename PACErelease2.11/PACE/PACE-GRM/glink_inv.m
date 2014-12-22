% inverse of the link function
function ginv=glink_inv(gfun,family,param)

    if strcmp(family,'binomial') == 1
        gfun(gfun<=0) = exp(-5);
        gfun(gfun>=1) = 1-exp(-5);
        ginv = log(gfun./(1-gfun));
    elseif strcmp(family,'poisson') == 1
        gfun(gfun<=0) = exp(-5);
        ginv = log(gfun);
    end
