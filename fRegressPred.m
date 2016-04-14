function ypredfd = fRegressPred(xfdcell, betaestcell, ybasisobj)
    % obtain predicted fd and matrix
    % given ybasisobj if the response variable is functional 
    % otherwise do not give input of ybasisobj, yredfd is a vector of scala
    
    % xfdcell ....... new functional covariates
    % betaestcell ... estimated beta from fRegress

    p = length(xfdcell);
    N = size(getcoef(xfdcell{1}), 2);
    
    if nargin == 3
    ynbasis = getnbasis(ybasisobj);
    rangeval = getbasisrange(ybasisobj);
    
    nfine = max(501, 10 * ynbasis + 1);
    tfine = linspace(rangeval(1), rangeval(2), nfine)';
    ymat = zeros(nfine, N);
    
    for j = 1:p
        xmat = eval_fd(tfine, xfdcell{j});
        betafdj = getfd(betaestcell{j});
        betavec = eval_fd(tfine, betafdj);
        ymat = ymat + xmat.*(betavec * ones(1, N));
    end
    ypredfd = smooth_basis(tfine, ymat, ybasisobj);
    
    elseif nargin == 2
        ymat = zeros(N, 1);
        for j=1:p
            xfdj = xfdcell{j};
            if isa_fd(xfdj)    
                xbasis  = getbasis(xfdj);
                xnbasis = getnbasis(xbasis);
                xrng    = getbasisrange(xbasis);
                nfine   = max(501,10*xnbasis+1);
                tfine   = linspace(xrng(1), xrng(2), nfine)';
                deltat  = tfine(2)-tfine(1);
                xmat    = eval_fd(tfine, xfdj);
                betafdj = getfd(betaestcell{j});
                betamat = eval_fd(tfine, betafdj);
                ymat = ymat + deltat.*(xmat'*betamat    - ...
                          0.5.*(xmat(1,    :)'*betamat(1)     + ...
                                xmat(nfine,:)'*betamat(nfine)));
            else
                betaj   = getcoef(getfd(betaestcell{j}));
                ymat = ymat + xfdj*betaj';
            end
        end
       ypredfd = ymat; 
       
    else
        error('Require either 2 or 3 input argument.');
    end
            
end