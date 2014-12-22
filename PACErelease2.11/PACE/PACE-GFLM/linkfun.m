% function for calculating g(eta),g^-1(mu),g'(eta) when g (inver of link
% function) is specified

function gv = linkfun(x,linktyp,ind);

if ind==0
    switch linktyp
        case 1
            gv=x;
        case 2
            gv=x.^2;
        case 3
            gv=log(x);
        case 4
            gv=log(x./(1-x));
        case 5
            gv=log(-log(1-x));
        case 6
            gv=1./x;
        case 7
            gv=sqrt(x);
    end
elseif ind==1
    switch linktyp
        case 1
            gv=ones(1,length(x));
        case 2
            gv=1./(2*x);
        case 3
            gv=x;
        case 4
            gv=x.*(1-x);
        case 5
            gv=(1-x).*(-log(1-x));
        case 6
            gv=-x.^2;
        case 7
            gv=2*sqrt(x);
    end
elseif ind==-1
    switch linktyp
        case 1
            gv=x;
        case 2
            gv=sqrt(x);
        case 3
            gv=exp(x);
        case 4
            gv=exp(x)./(1+exp(x));
        case 5
            gv=1-exp(-exp(x));
        case 6
            gv=1./x;
        case 7
            gv=x.^2;
    end   
end

end