function  ss = inprod_TPbasis2(basis1,  basis2,  basis3,  basis4, ...
                              Lfdobj1, Lfdobj2, Lfdobj3, Lfdobj4, ...
                              rng, wtfd, EPS, JMAX, JMIN)
%  INPROD_TPBASIS  Computes vectorized four-tensor of inner products of the
%    tensor product of values of respective linear differential operators
%    applied to four bases.  
%    The inner products are approximated by numerical integration using  
%    Romberg integration with the trapezoidal rule. 
%    The inner s can be over a reduced range in RNG and use a 
%    scalar weight function in WTFD.
%
%  Arguments:
%  BASIS1, BASIS2, BASIS3 and BASIS4 ...  these are basis objects,
%            and the inner products of all quadruples of 
%            BASIS1, BASIS2, BASIS3 and BASIS4 are computed.
%  LFDOBJ1, LFDOBJ2, LFDOBJ3 and LFDOBJ4 ... linear differential operators  
%                    for the corresponding basis functions.
%  RNG  ...  Limits of integration
%  WTFD ...  A functional data object defining a weight
%  EPS  ...  A convergence criterion, defaults to 1e-4.
%  JMAX ...  Maximum number of Richardson extrapolation iterations.
%            Defaults to 16.
%  JMIN ...  Minimum number of Richardson extrapolation iterations.
%            Defaults to 5.
%
%  Return:
%  An order NBASIS1*NBASIS2*NBASIS3*NBASIS4 matrix SS of inner products 
%  for each possible pair quadruple of basis functions.

%  last modified 6 February 2015

%  set up default values of arguments

if nargin < 13 || isempty(JMIN),    JMIN = 5;             end
if nargin < 12 || isempty(JMAX),    JMAX = 16;            end
if nargin < 11 || isempty(EPS),     EPS  = 1E-6;          end
if nargin < 10 || isempty(wtfd),    wtfd = 0;             end
if nargin <  9 || isempty(rng),     rng  = [];            end
if nargin <  8 || isempty(Lfdobj4), Lfdobj4 = int2Lfd(0); end
if nargin <  7 || isempty(Lfdobj3), Lfdobj3 = int2Lfd(0); end
if nargin <  6 || isempty(Lfdobj2), Lfdobj2 = int2Lfd(0); end
if nargin <  5 || isempty(Lfdobj1), Lfdobj1 = int2Lfd(0); end

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

ncum  = cumprod([nbasis1, nbasis2, nbasis3, nbasis4]);
nprod = ncum(4);

range1  = getbasisrange(basis1);
range2  = getbasisrange(basis2);
range3  = getbasisrange(basis3);
range4  = getbasisrange(basis4);

if nargin < 9 || isempty(rng), rng = range1; end
if rng(1) < range1(1) || rng(2) > range1(2) || ...
   rng(1) < range2(1) || rng(2) > range2(2) || ...
   rng(1) < range3(1) || rng(2) > range3(2) || ...
   rng(1) < range4(1) || rng(2) > range4(2)
      error('Limits of integration are inadmissible.');
end

%  set up first iteration using only boundary values

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h     = ones(JMAXP,1);
h(2)  = 0.25;
s     = zeros(JMAXP,nprod);
if ~isnumeric(wtfd)
    wtvec = eval_fd(wtfd, rng);
else
    wtvec = ones(2,1);
end
bmat1  = eval_basis(rng, basis1, Lfdobj1);
bmat2  = eval_basis(rng, basis2, Lfdobj2);
bmat3  = eval_basis(rng, basis3, Lfdobj3);
bmat4  = eval_basis(rng, basis4, Lfdobj4);
tensorprod = zeros(nprod,1);
for m=1:2
    for i=1:nbasis1
        if bmat1(m,i) ~= 0
            for j=1:nbasis2
                if bmat2(m,j) ~= 0
                    for k=1:nbasis3
                        if bmat3(m,k) ~= 0
                            for l=1:nbasis4
                                if bmat4(m,l) ~= 0
                                    ijkl = (l-1)*ncum(3) + ...
                                           (k-1)*ncum(2) + ...
                                           (j-1)*ncum(1) + i;
                                    tensorprod(ijkl) = ...
                                        tensorprod(ijkl) + ...
                                        bmat1(m,i)*bmat2(m,j)* ...
                                        bmat3(m,k)*bmat4(m,l)*wtvec(m);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
chs    = width.*tensorprod'./2;
s(1,:) = chs;
tnm    = 0.5;

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
    tensorprod = zeros(nprod,1);
    for m=1:nx
        for i=1:nbasis1  
            if bmat1(m,i) ~= 0
                for j=1:nbasis2  
                    if bmat2(m,j) ~= 0
                        for k=1:nbasis3 
                            if bmat3(m,k) ~= 0
                                for l=1:nbasis4  
                                    if bmat4(m,l) ~= 0
                                        ijkl = (l-1)*ncum(3) + ...
                                               (k-1)*ncum(2) + ...
                                               (j-1)*ncum(1) + i;
                                        tensorprod(ijkl) = ...
                                            tensorprod(ijkl) + ...
                                            bmat1(m,i)*bmat2(m,j)* ...
                                            bmat3(m,k)*bmat4(m,l)*wtvec(m);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    chs = width.*tensorprod'./tnm;
    chsold = s(iter-1,:);
    s(iter,:) = (chsold + chs)./2;
    if iter >= 5
        ind   = (iter-4):iter;
        ya    = s(ind,:);
        xa    = h(ind);
        absxa = abs(xa);
        [absxamin, ns] = min(absxa);
        cs = ya;
        ds = ya;
        y  = ya(ns,:);
        ns = ns - 1;
        for m = 1:4
            for i = 1:(5-m)
                ho      = xa(i);
                hp      = xa(i+m);
                w       = (cs(i+1,:) - ds(i,:))./(ho - hp);
                ds(i,:) = hp.*w;
                cs(i,:) = ho.*w;
            end
            if 2*ns < 5-m
                dy = cs(ns+1,:);
            else
                dy = ds(ns,:);
                ns = ns - 1;
            end
            y = y + dy;
        end
        ss = y;
        errval = max(abs(dy));
        ssqval = max(abs(ss));
        if all(ssqval > 10*eps)
            crit = errval./ssqval;
        else
            crit = errval;
        end
        % disp(['Iteration ',num2str(iter),'  Crit = ',num2str(crit)])
        if crit < EPS && iter >= JMIN
            ss = ss';
            % disp(['Iteration ',num2str(iter)])
            % disp(['Crit = ',num2str(crit)])
            return
        end
    end
    s(iter+1,:) = s(iter,:);
    h(iter+1)   = 0.25.*h(iter);
end
disp(['No convergence after ',num2str(JMAX),' steps in INPROD_TPBASIS.']);

