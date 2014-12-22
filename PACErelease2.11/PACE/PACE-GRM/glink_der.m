% derivative of the link function
function gder = glink_der(x,family,param)

    if strcmp(family,'binomial') == 1
        gder = exp(x)./(1+exp(x)).^2;
    elseif strcmp(family,'poisson') == 1
        gder = exp(x);
    end
