function  ss = inprod_TPbasis(basis1,  basis2,  basis3,  basis4, ...
                              Lfdobj1, Lfdobj2, Lfdobj3, Lfdobj4, ...
                              rng, wtfd, EPS, JMAX, JMIN)
%  INPROD_TPBASIS  Computes matrix of inner products of the
%    tensor product of four bases, by numerical integration using  
%       Romberg integration with the trapezoidal rule.
%
%  Arguments:
%  BASIS1, BASIS2, BASIS3 and BASIS4 ...  these are basis objects,
%            and the inner products of (BASIS1, BASIS2) and (BASIS3,BASIS4)
%            are computed.
%  LFDOBJ1, LFDOBJ2, LFDOBJ13 and Lfdobj24 ...  differential operators for 
%                    respective basis funcitons
%  RNG  ...  Limits of integration
%  WTFD ...  A functional data object defining a weight
%  EPS  ...  A convergence criterion, defaults to 1e-4.
%  JMAX ...  Maximum number of Richardson extrapolation iterations.
%            Defaults to 15.
%  JMIN ...  Minimum number of Richardson extrapolation iterations.
%            Defaults to 5.
%
%  Return:
%  An order NREP1*NREP2 matrix SS of inner products for each possible pair
%  of basis functions.

%  last modified 30 June 2014

%  set up default values of arguments

if nargin < 13, JMIN = 5;             end
if nargin < 12, JMAX = 15;            end
if nargin < 11, EPS  = 1E-4;          end
if nargin < 10, wtfd = 0;             end
if nargin <  9, rng  = [];            end
if nargin <  8, Lfdobj4 = int2Lfd(0); end
if nargin <  7, Lfdobj3 = int2Lfd(0); end
if nargin <  6, Lfdobj2 = int2Lfd(0); end
if nargin <  5, Lfdobj1 = int2Lfd(0); end

%  check LFD objects

Lfdobj1 = int2Lfd(Lfdobj1);
Lfdobj2 = int2Lfd(Lfdobj2);
Lfdobj3 = int2Lfd(Lfdobj3);
Lfdobj4 = int2Lfd(Lfdobj4);
  
%  check WTFD

if  isa_fd(wtfd)
    coefw = getcoef(wtfd);
    coefd = size(coefw);
    if coefd(2) > 1
        error('Argument WTFD is not a single function');
    end
end
 
%  check basis function objects

if ~(isa_basis(basis1) && isa_basis(basis2) && ...
     isa_basis(basis3) && isa_basis(basis4))
    error ('The four first arguments are not basis objects.');
end

%  determine NBASIS1 and NASIS2, and check for common range

nbasis1 = getnbasis(basis1) - length(getdropind(basis1));
nbasis2 = getnbasis(basis2) - length(getdropind(basis2));
nbasis3 = getnbasis(basis3) - length(getdropind(basis3));
nbasis4 = getnbasis(basis4) - length(getdropind(basis4));
range1  = getbasisrange(basis1);
range2  = getbasisrange(basis2);
range3  = getbasisrange(basis3);
range4  = getbasisrange(basis4);
if nargin < 9 || isempty(rng), rng = range1; end
if rng(1) < range1(1) || rng(2) > range1(2) || ...
   rng(1) < range2(1) || rng(2) > range2(2) || ...
   rng(1) < range3(1) || rng(2) > range3(2) || ...
   rng(1) < range4(1) || rng(2) > range4(2)
    disp(rng)
    disp(range1)
    disp(range2)
    disp(range3)
    disp(range4)
      error('Limits of integration are inadmissible.');
end

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h     = ones(JMAXP,1);
h(2)  = 0.25;
s = reshape(zeros(JMAXP*nbasis1*nbasis2*nbasis3*nbasis4,1),...
            [JMAXP,nbasis1,nbasis2,nbasis3,nbasis4]);
%  the first iteration uses just the endpoints
if ~isnumeric(wtfd)
    wtvec = eval_fd(wtfd, rng);
else
    wtvec = ones(2,1);
end
bmat1  = eval_basis(rng, basis1, Lfdobj1);
bmat2  = eval_basis(rng, basis2, Lfdobj2);
bmat3  = eval_basis(rng, basis3, Lfdobj3);
bmat4  = eval_basis(rng, basis4, Lfdobj4);
tensorprod = zeros(nbasis1,nbasis2,nbasis3,nbasis4);
for i=1:nbasis1
    for j=1:nbasis2
        for k=1:nbasis3
            for l=1:nbasis4
                for m=1:2
                    tensorprod(i,j,k,l) = tensorprod(i,j,k,l) + ...
                                        bmat1(m,i)*bmat2(m,j)* ...
                                        bmat3(m,k)*bmat4(m,l)*wtvec(m);
                end
            end
        end
    end
end
chs = width.*tensorprod./2;
s(1,:,:,:,:) = chs;
tnm = 0.5;

%  now iterate to convergence

for iter = 2:JMAX
    % disp(iter)
    tnm = tnm.*2;
    del = width./tnm;
    x   = rng(1)+del/2:del:rng(2);
    nx  = length(x);
    if ~isnumeric(wtfd)
        wtvec = eval_fd(wtfd, x);
    else
        wtvec = ones(nx,1);
    end
    bmat1  = eval_basis(x, basis1, Lfdobj1);
    bmat2  = eval_basis(x, basis2, Lfdobj2);
    bmat3  = eval_basis(x, basis3, Lfdobj3);
    bmat4  = eval_basis(x, basis4, Lfdobj4);
    tensorprod = zeros(nbasis1,nbasis2,nbasis3,nbasis4);
    for i=1:nbasis1
        for j=1:nbasis2
            for k=1:nbasis3
                for l=1:nbasis4
                    for m=1:nx
                        tensorprod(i,j,k,l) = tensorprod(i,j,k,l) + ...
                                            bmat1(m,i)*bmat2(m,j)* ...
                                            bmat3(m,k)*bmat4(m,l)*wtvec(m);
                    end
                end
            end
        end
    end
    chs = width.*tensorprod./tnm;
    chsold = reshape(s(iter-1,:,:,:,:),size(chs));
    s(iter,:,:,:,:) = (chsold + chs)./2;
    if iter >= 5
        ind = (iter-4):iter;
        ya = s(ind,:,:,:,:);
        xa = h(ind);
        absxa = abs(xa);
        [absxamin, ns] = min(absxa);
        cs = ya;
        ds = ya;
        y  = squeeze(ya(ns,:,:,:,:));
        ns = ns - 1;
        for m = 1:4
            for i = 1:(5-m)
                ho      = xa(i);
                hp      = xa(i+m);
                w       = (cs(i+1,:,:,:,:) - ds(i,:,:,:,:))./(ho - hp);
                ds(i,:,:,:,:) = hp.*w;
                cs(i,:,:,:,:) = ho.*w;
            end
            if 2*ns < 5-m
                dy = squeeze(cs(ns+1,:,:,:,:));
            else
                dy = squeeze(ds(ns,:,:,:,:));
                ns = ns - 1;
            end
            y = y + dy;
        end
        ss = reshape(y, nbasis1,nbasis2,nbasis3,nbasis4);
        errval = max(abs(reshape(dy,nbasis1*nbasis2*nbasis3*nbasis4,1)));
        ssqval = max(abs(reshape(ss,nbasis1*nbasis2*nbasis3*nbasis4,1)));
        if all(ssqval > 0)
            crit = errval./ssqval;
        else
            crit = errval;
        end
        if crit < EPS && iter >= JMIN
            ss = squeeze(ss);
            return
        end
    end
    s(iter+1,:,:,:,:) = s(iter,:,:,:,:);
    h(iter+1)   = 0.25.*h(iter);
end
disp(['No convergence after ',num2str(JMAX),' steps in INPROD_BASIS.']);

