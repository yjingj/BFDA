function [Bmat, mu, Fval] = ...
             glm_mda(yCell, Zmat, distCell, lamRmat, Wtvec, Bmat0, addterm)
%GLM_MDA Fits a generalized linear model with regularization to a
%  N by M matrix of data, each column having a potentially different
%  GLM model.
%  Arguments:
%
%  YCELL    ... May be:
%                 an N by P matrix of data to be fitted or a cell array
%                 of length P.  
%                 For the matrix format, the binomial data case for any
%                 column is not permitted, but a column can contain only
%                 0's and 1's, in which case the glm family is 'binary'.
%                 If one or more of the observations are in the binomial
%                 family, then the celld array format is required because
%                 in the binomial case with local sample sizes M_i,
%                 a cell array of length 2 is required in the corresponding
%                 cell of YCELL, the first cell of which contains the 
%                 observed frequencies, and the second of which contains  
%                 the corresponding sample sizes. 
%  ZMAT     ... An N by Q matrix of covariates or independent variable
%               values.
%  DISTCELL ... A cell array of length P, each cell of which contains
%               a string indicating which of the five GLM family members
%               is assumed:
%                 'normal' or 'gaussian' or 'Gaussian'
%                 'binomial' or 'binary' or 'Bernoulli'
%                 'poisson'
%                 'gamma'
%                 'inverse gaussian' or 'inverse Gaussian'
%  LAMRMAT  ... a \lambda R, that is, a roughness penalty matrix R of 
%               order equal to the number of basis functions used or number
%               of columns of Zmat multiplied by a scalar roughness 
%               penalty parameter \lambda
%  WTVEC    ... a vector of prior weights, such as the inverses of the
%               relative variance of each observation.
%  Bmat0    ... starting values for regression coefficients
%  ADDTERM  ... a addterm with a coefficient fixed at 1.0.
%
%  Returns:
%  BMAT     ... Final estimate of regression coefficients
%  MU       ... Matrix of estimated expectations
%  FVAL     ... Penalized deviance value
%
%   Last modified 18 November 2014 by Jim Ramsay

%--------------------------------------------------------------------------
%                    Check arguments
%--------------------------------------------------------------------------

if nargin < 3
    error('Number of arguments is less than 3.');
end

%  set default argument values

if nargin < 5, Wtvec    = []; end
if nargin < 6, Bmat0    = []; end
if nargin < 7, addterm  = []; end

if nargin < 4
    error('Less than four arguments for function glm_fda.');
end
   
%  check dimensions of ZMAT and YCELL

[N, q] = size(Zmat);
if iscell(yCell)
    %  YCELL is a cell array of length P
    p = length(yCell);
    %  check dimension of object in each cell
    Ymat = zeros(N,p);
    for j=1:p
        yCellj = yCell{j};
        if iscell(yCellj)
            %  object is a cell array, which must of length 2
            %  for binomial observations
            if length(yCellj) ~= 2
                error(['A cell of YCELL contains a cell array which', ...
                       ' is not of length 2.']);
            end
            Ntmpj = size(yCellj{1},1);
            if N ~= Ntmpj
                error(['ZMAT and a cell of YCELL do not have', ...
                    ' the same number of rows.']);
            end
            Ymat(:,j) = yCellj{1};
        else
            Ntmpj = size(yCellj,1);
            Ymat(:,j) = yCellj;
        end
        if N ~= Ntmpj
            error(['ZMAT and a cell of YCELL do not have', ...
                ' the same number of rows.']);
        end
    end
else
    [Ntmp, p] = size(yCell);
    if N ~= Ntmp
        error(['ZMAT and matrix YCELL do not have', ...
               ' the same number of rows.']);
    end
    %  convert matrix to a cell array for code simplicity
    temp = yCell;
    yCell = cell(1,p);
    for j=1:p
        yCell{j} = temp(:,j);
    end
end

%  define default weight vector WTVEC and check for positivity

if isempty(Wtvec)
    Wtvec = ones(N,1);
end

if any(Wtvec <= 0)
    error('Non-positive values of WTVEC found.');
end

%--------------------------------------------------------------------------
%  Process DISTCELL and define anonymous functions according to the 
%  distribution of each cell in YCELL.  A struct object is set up for each.
%     devFn   ... the deviance or loss function, 
%                 called after convergence is achieved
%     stdFn   ... the scale factor multiplying D eta
%                 called second inside loop        
%     linkFn  ... link function, eta = linkFn(mu),
%                 called prior to loop, maps data space into real line
%     DlinkFn ... derivative of the link function wrt to mu
%                 called first inside loop
%     IlinkFn ... the inverse of the link function IlinkFn[eta] = mu,
%                 called last inside loop, maps eta into data space
% Then set a starting value for the mean mu, avoiding boundary values.
%--------------------------------------------------------------------------

if ~iscell(distCell)
    error('DISTCELL is not a cell array.');
end
if length(distCell) ~= p
    error('DISTCELL is not same length as YCELL.');
end

%  Process each cell of DISTCELL, and in each set up a struct object
%  containing required anonymous functions and initial value for MU
for j=1:p
    diststrj = distCell{j};
    switch diststrj
        case 'normal'
            %  Note:  Ymat can be any real number, no restrictions
            Yvec    = yCell{j};
            devFn   = @(mu,Yvec) (Yvec - mu).^2;
            stdFn   = @(mu)  ones(size(mu));
            linkFn  = @(mu)  mu;
            DlinkFn = @(mu)  ones(size(mu));
            IlinkFn = @(eta) eta;
            mu      = Yvec;
        case 'binomial'
            if isnumeric(Yvec)
                %  If Yvec is binary, M is taken to be 1 (set below)
                %  and it must be a binary matrix containing only 
                %  0's and 1's
                if any(Yvec ~= 0 & Yvec ~= 1)
                    error(['For binomial case, Yvec a single column ', ...
                        ' but contains values other than 0 or 1.']);
                end
                M = ones(N,p);
            elseif iscell(Yvec) && length(Yvec) == 2
                %  If Yvec is a cell array of length 2, then first cell 
                %  contains a matrix containing the number of successes and 
                %  the second cell either contains a matrix of the same 
                %  size as the matrix in Yvec{1} or a single positive 
                %  integer.  
                %  These values or this value is the number of trials M
                %  for a binomial or bernoulli distribution.
                %  M must be a positive integer.
                Freq = Yvec{1};
                M    = Yvec{2};
                if length(M) == 1
                    M = M*ones(N,p);
                end
                if ~all(size(M) == size(Freq))
                    error(['DISTR is binomial and matrix M is not the', ...
                        ' same size as matrix FREQ']);
                end
                if any(any(M < 0))
                    error(['DISTR is binomial and one or more values', ...
                        ' in M have nonpositive values']);
                end
                if any(any(floor(M) ~= M))
                    error(['DISTR is binomial and one or more values', ...
                        ' in M have noninteger values.']);
                end
                %  Redefine Yvec is the proportion of sucesses
                Yvec = Freq./M;
            else
                error(['YVEC has incorrect dimensions for binomial', ...
                    ' distributin, or is of wrong type.']);
            end
            devFn   = @(mu,Yvec) ...
                2*M.*(Yvec.*log((Yvec+(Yvec==0))./mu) + ...
                (1-Yvec).*log((1-Yvec+(Yvec==1))./(1-mu)));
            stdFn   = @(mu)  sqrt(mu.*(1-mu)./M);
            linkFn  = @(mu)  log(mu./(1-mu));
            DlinkFn = @(mu)  1./(mu.*(1-mu));
            loBnd   = log(eps);
            upBnd   = -loBnd;
            IlinkFn = @(eta) 1./(1 + exp(-constrain(eta,loBnd,upBnd)));
            mu      = (M.*Yvec + 0.5)./(M + 1);
        case 'poisson'
            %  Note: Yvec must not contain negative numbers
            if any(any(Yvec < 0))
                error(['DISTR is poisson and Yvec contains negative ', ...
                       'values']);
            end
            devFn   = @(mu,Yvec) 2*(Yvec.*(log((Yvec+(Yvec==0))./mu)) - ...
                (Yvec - mu));
            stdFn   = @(mu)  sqrt(mu);
            linkFn  = @(mu)  log(mu);
            DlinkFn = @(mu)  1./mu;
            loBnd   = log(eps);
            upBnd   = -loBnd;
            IlinkFn = @(eta) exp(constrain(eta,loBnd,upBnd));
            mu      = Yvec + 0.25;
        case 'gamma'
            %  Note:  Yvec must contain only positive numbers
            if any(Yvec <= 0)
                error('DISTR is gamma and Y contains nonpositive values');
            end
            devFn   = @(mu,Yvec) 2*(-log(Yvec./mu) + (Yvec - mu)./mu);
            stdFn   = @(mu) mu;
            linkFn  = @(mu)  1./mu;
            DlinkFn = @(mu) -1./mu.^2;
            loBnd   = eps;
            upBnd   = 1/loBnd;
            IlinkFn = @(eta) 1./constrain(eta,loBnd,upBnd);
            mu      = max(Yvec, eps);
        case 'inverse gaussian'
            %  Note:  Yvec must contain only positive numbers
            if any(any(Yvec <= 0))
                error(['DISTR is inverse gaussian and Y contains ', ...
                    'nonpositive values']);
            end
            devFn   = @(mu,Yvec) ((Yvec - mu)./mu).^2./ Yvec;
            stdFn   = @(mu)  mu.^(3/2);
            loBnd   = eps.^(1/2);
            upBnd   = 1/loBnd;
            linkFn  = @(mu)  constrain(mu,loBnd,upBnd).^(-2);
            DlinkFn = @(mu)  -2*mu.^(-3);
            IlinkFn = @(eta) constrain(eta,loBnd,upBnd).^(-1/2);
            mu      = Yvec;
        otherwise
            error('Distribution name is invalid.');
    end
    distStructj.family  = diststrj;
    distStructj.devFn   = devFn;
    distStructj.stdFn   = stdFn;
    distStructj.loBnd   = lobnd;
    distStructj.upBnd   = upBnd;
    distStructj.linkFn  = linkFn;
    distStructj.DlinkFn = DlinkFn;
    distStructj.IlinkFn = IlinkFn;
    distStructj.mu      = mu;
    distCell{j} = distStructj;
end

%--------------------------------------------------------------------------
%                   Initialize mu and eta from Yvec.
%--------------------------------------------------------------------------

% compute eta = E(y) from mu

eta  = zeros(N,p);

for j=1:p
    distStructj = distCell{j};
    for i=1:N
        eta(i,j) = distStructj.linkFni(distStructj.mu(i));
    end
end

%--------------------------------------------------------------------------
%                        Set up for iterations
%--------------------------------------------------------------------------

iter     = 0;
iterLim  = 100;
convcrit = 1e-6;
sqrtwt   = sqrt(Wtvec);

%  set up starting value Bmat0 if required

if isempty(Bmat0)
    Bmat0 = zeros(q,p);
end
Bmat = Bmat0;

% Enforce limits on mu to guard against an inverse linkFn that doesn't map 
% into the support of the distribution.

for j=1:p
    distStructj = distCell{j};
    switch distStructj.familu
        case 'binomial'
            % mu is a probability, so order one is the natural scale, 
            % and eps is areasonable lower limit on that scale 
            % (plus it's symmetric).
            muLims = [eps 1-eps];
        case {'poisson' 'gamma' 'inverse gaussian'}
            % Here we don't know the natural scale for mu, so make the 
            % lower limit small.  This choice keeps mu^4 from underflowing.  
            %  No upper limit.
            muLims = realmin.^.25;
    end
end

%--------------------------------------------------------------------------
%                       Start of GLM iteration loop
%--------------------------------------------------------------------------

while iter <= iterLim
    iter = iter+1;

    % Compute adjusted dependent variable for least squares fit
    
    Deta = zeros(N,p);
    stdm = zeros(N,p);
    mu   = zeros(N,p);
    for j=1:p
        distStructj = distCell{j};
        for i=1:N
            mu(i,j)   = distStructj.mu(i);
            Deta(i,j) = distStructj.DlinkFni(mui);
            stdm(i,j) = distStructj.stdFni(mui);
        end
    end
    Umat = eta + (Ymat - mu).*Deta;
  
    % Compute IRLS weights the inverse of the variance function
    
    sqrtw = (sqrtwt*ones(1,p))./(abs(Deta).*stdm);
    
    % Compute coefficient estimates for this iteration - the IRLS step
    
    Bmat_old   = Bmat;
    if ~isempty(addterm)
        Ymattmp = Umat - addterm;
    else
        Ymattmp = Umat;
    end
    Ymatw   = Ymattmp.*sqrtw;
    Zmatw   = Zmat.*(sqrtwt*ones(1,nbasis));
    if isempty(lamRmat)
        Mmat = Zmatw'*Zmatw;
    else
        Mmat = Zmatw'*Zmatw + lamRmat;
    end
    Bmat    = Mmat\Zmatw'*Ymatw;
    if ~isempty(addterm)
        eta = Zmat*Bmat + addterm;
    else
        eta = Zmat*Bmat;
    end
    for j=1:p
        distStructj = distCell{j};
        for i=1:N
            mu(i,j) = distStructj.IlinkFni(eta(i,j));
        end
    end

    % Force mean in bounds, in case the linkFn function is faulty
    
    for j=1:p
        distStructj = distCell{j};
        for i=1:N
            switch distStructj.family
                case 'binomial'
                    if mu(i,j) < muLims(1) || muLims(2) < mu(i,j)
                        mu(i,j) = max(min(mu(i,j),muLims(2)),muLims(1));
                    end
                case {'poisson' 'gamma' 'inverse gaussian'}
                    if mu(i,j) < muLims(1)
                            mu(i,j) = max(mu(i,j),muLims(1));
                    end
            end
        end
    end

    % Check stopping conditions
   
    if max(max(abs(Bmat-Bmat_old))) < ...
            convcrit*max(max(abs(Bmat_old))) 
        break; 
    end
    
end

%--------------------------------------------------------------------------
%                    end of GLM iteration loop
%--------------------------------------------------------------------------

if iter > iterLim
    warning(['smooth_GLM:','Iteration'],'Iteration limit reached.');
end

if nargout > 3
    % Sum components of deviance to get the total deviance.
    Deviance = zeros(N,p);
    for j=1:p
        distStructj = distCell{j};
        for i=1:N;
            di = distStructj.devFn(mu(i,j),Ymat(i,j));
            Deviance(i,j) = Wtvec(i)*di;
        end
    end
end

Fval = sum(sum(Deviance));
if ~isempty(lamRmat)
    Fval = Fval + trace(Bmat'*lamRmat*Bmat);
end

