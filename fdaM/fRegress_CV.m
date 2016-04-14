function [SSE_CV, errfd] = ...
              fRegress_CV(yfdPar, xfdcell, betacell, wt, CVobs)
% FREGRESS_CV computes cross-validated error sum of squares
% for scalar or functional responses. NOTE: ordinary and
% generalized cross validation scores are now returned by fRegress
% when scalar responses are used.

%  Last modified 28 July 2012 by Jim Ramsay

if nargin < 4, wt = [];  end

%  check the arguments

[yfdPar, xfdcell, betacell, wt] = ...
               fRegress_argcheck(yfdPar, xfdcell, betacell, wt);
  
%  extract dimensions of the data and the analysis

p = length(xfdcell);
N = size(getcoef(xfdcell{1}),2);
%  default value for CVobs
if nargin < 5 || isempty(CVobs), CVobs = 1:N;  end
M = length(CVobs);

%  branch to either scalar or functional dependent variable

if isnumeric(yfdPar) 
    
    %  Dependent variable is scalar
    
    yvec   = yfdPar;
    SSE_CV = 0;
    errvec = zeros(M,1);
    for m=1:M
        i        = CVobs(m);
        indexi   = find((1:N) ~= i);
        if ~isempty(wt)
            wti  = wt(indexi);
        else
            wti  = [];
        end
        xfdcelli = cell(p,1);
        for j=1:p
            xfd = xfdcell{j};
            xfdcelli{j} = xfd(indexi);
        end
        yveci         = yvec(indexi);
        fRegressStri = fRegress(yveci, xfdcelli, betacell, wti);
        betaestcelli  = fRegressStri.betahat;
        yhati = 0;
        for j=1:p
            betafdj = getfd(betaestcelli{j});
            xfdj    = xfdcell{j};
            bbasisj = getbasis(betafdj);
            rangej  = getbasisrange(bbasisj);
            nfine   = max(501, getnbasis(bbasisj)*10+1);
            tfine   = linspace(rangej(1), rangej(2), nfine)';
            delta   = tfine(2)-tfine(1);
            betavec = eval_fd(tfine, betafdj);
            xveci   = eval_fd(tfine, xfdj(i));
            yhati   = yhati + delta.*(sum(xveci.*betavec) - ...
                0.5.*( xveci(1)    *betavec(1) + ...
                       xveci(nfine)*betavec(nfine) ));
        end
        errvec(i) = yvec(i) - yhati;
        SSE_CV    = SSE_CV + errvec(i).^2;
    end
    errfd = errvec;
else
    
    %  Dependent variable is functional
    
    yfd      = getfd(yfdPar);
    SSE_CV   = 0;
    errcoefs = [];
    for m = 1:M
        %  index of case to eliminate and indices of remaining
        i = CVobs(m);
        indexi = find((1:N) ~= i);
        %  eliminate case i from the weights
        if ~isempty(wt)
            wti  = wt(indexi);
        else
            wti  = [];
        end
        %  eliminate case i from covariates
        xfdcelli = cell(p,1);
        for j=1:p
            xfd = xfdcell{j};
            xfdcelli{j} = xfd(indexi);
        end
        yfdi = yfd(indexi); 
        %  carry out the functional regression analysis
        fRegresscelli = fRegress(yfdi,xfdcelli,betacell,wti);
% fRegressStruct = ...
%         fRegress(yfdPar, xfdcell, betacell, wt, y2cMap, SigmaE)
% yfdPar=yfdi;
% xfdcell=xfdcelli;
% wt=wti;
        %  extract the regression coefficient functions
        betaestcelli = fRegresscelli.betahat;
        %  compute the fit to the data for case i
        yhatfdi = 0;
        for j=1:p
            betafdParj = betaestcelli{j};
            betafdj    = getfd(betafdParj);
            xfdj       = xfdcell{j};
            xfdij      = xfdj(i);
            tempfd     = xfdij.*betafdj;
            yhatfdi    = yhatfdi + tempfd;
        end
        %  compute the residual function
        errfdi   = yfd(i) - yhatfdi;
        %  increment the error sum of squares by the integral of the
        %  square of the residual function
        SSE_CV   = SSE_CV + inprod(errfdi,errfdi);
        %  add the coefficients for the resiudal function
        errcoefs = [errcoefs, getcoef(errfdi)];
    end
    %  set up the functional data object for the residual fns
    if nargout > 1
        errfd = fd(errcoefs,getbasis(errfdi));
    end
end



