% function beta = transferX(xx, nder, gridin, gridout)
% Calculates transfer functions of Gaussian processes X. 
% \beta(s,t) = E(X^{(\nu)}(t)|X(s))
% Reference: M\"uller, H.G. and Yang, W. (2009). Dynamic relations for
% sparsely sampled Gaussian processes. TEST.
%======
%Input:
%======
%  xx:      the returned values from FPCder. See FPCder() for more details.
%  nder:    the order \nu derivative of X at t that is regressed on X(s).
%  gridin:  optional, the input grid of the transfer function, corresponds to argument s.
%  gridout: optional, the output grid of the transfer function, corresponds to argument t.
%
%=======
%Output:
%=======
%  beta:    transfer function, a matrix whose dimension is determined by
%           gridin and gridout.
function beta = transferX(xx, nder, gridin, gridout)
    if nargin < 4
        if nargin == 2
            gridin = getVal(xx,'out21copy');
        else
            fprintf(1,'Error: Not enough input arguments!');
            beta = [];
            return;
        end
        gridout = gridin;
    end
    lambda = getVal(xx,'lambda');
    phi = getVal(xx,'phicopy');
    if ~iscell(phi) && nder == 0
        phider = phi;
    elseif ~iscell(phi) && nder ~= 0
        fprintf(1,'Error: xx must be the output from FPCder()!');
        beta = [];
        return;
    else
        phider = phi{nder+1};
        phi = phi{1};
    end
    out1 = getVal(xx,'out1copy');
    k = size(lambda,2);
    phis = interp1(out1,phi, gridin,'spline');
    phit = interp1(out1,phider, gridout,'spline');
    beta = zeros(length(gridin),length(gridout));
    if size(phis,1) == 1
        phis = phis';
    end
    if size(phit,1) == 1
        phit = phit';
    end
    for i = 1:k
        beta = beta + lambda(i)*repmat(phis(:,i),1,length(gridout)).* repmat(phit(:,i)',length(gridin),1);
    end
    beta = beta./repmat((phis.^2)*lambda',1,length(gridout));
end
