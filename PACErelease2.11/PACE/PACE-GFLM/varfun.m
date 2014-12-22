% function for calculating variance when variance function is specified

function sigma2 = varfun(mu,Y,vartyp,theta);

switch vartyp
    case 1
        sigma2=repmat(var(Y-mu),[1,length(mu)]);
    case 2 
        sigma2=mu.*(1-mu)./theta;
    case 3
        sigma2=mu;
    case 4
        sigma2=mu.^2./theta;
end

end