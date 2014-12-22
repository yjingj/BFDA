% quasi-likelihood 
function ql = QL_param(x,y,family,param,sigma)

    if strcmp(family,'binomial') == 1
        ql = sum(log(x.^y)+log((1-x).^(1-y)))/sigma;
    elseif strcmp(family,'poisson') == 1
        ql = sum(log(x.^y)-x-log(factotial(y)))/sigma;
    end

    