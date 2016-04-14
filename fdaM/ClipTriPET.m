function [p,e,t] = ClipTriPET(T, M, B)
% ClipTriPET constructs a p e t mesh on the parallelogram
%          [(0,0), (T,0), (0,T), (-T,T)]
% The left boundary runs from (T-B,0) to (-B,T) and there are
%  M intervals on each side, with the left half of the parallelogram
%  being empty.
%
%  Last modified 7 June by Jim

%  check T

if T <= 0
    error('T is not positive.');
end

%  set default values for M and B

if nargin < 2, M = 1;  end
if nargin < 3, B = M;  end

%  check B

% if B < 1 || B > M
%     error('Incorrect value for B');
% end

%  compute points array

npts = (B+1)*(M+1);
ntri = 2*B*M;
p = zeros(2,npts);
X = linspace(-B*T/M,T,M+B+1);
m2p  = 0;
indi = 1:M+1;
for i=1:B+1
    m1p = m2p + 1;
    m2p = m2p + M + 1;
    p(1,m1p:m2p) = X(indi+B-i+1);
    p(2,m1p:m2p) = X(M+B+2-indi);
end

%  compute edge indices

nedg = 2*(M+B);
e = zeros(7,nedg);
e(6,:) = 1;
e(5,:) = 1:nedg;
ind1 = 1:M;
ind2 = 1:B;
% right boundary
e(1,ind1) = 1:M;
e(2,ind1) = 2:M+1;
% bottom boundary
e(1,ind2+M) = (M+1):(M+1):B*(M+1);
e(2,ind2+M) = 2*(M+1):(M+1):(B+1)*(M+1);
% left boundary
e(1,ind1 + B+M) =  npts:   (-1):(B*(M+1)+2);
e(2,ind1 + B+M) = (npts-1):(-1):(B*(M+1)+1);
% top boundary
e(1,ind2+B+2*M) = (B*(M+1)+1):-(M+1):(M+2);
e(2,ind2+B+2*M) = ((B-1)*(M+1)+1):-(M+1):1;

%  compute triangle indices

m2t = 0;
m2p = 0;
t = zeros(4,ntri);
t(4,:) = 1;
for i=1:B
    m1t = m2t +  1;
    m2t = m2t + 2*M;
    m1p = m2p + 1;
    m2p = m2p + M + 1;
    indup = m1p:m2p-1;
    inddn = M+1+(m1p+1:m2p);
    indtri1 = m1t:2:(m2t-1);
    t(1,indtri1) = inddn;
    t(2,indtri1) = indup;
    t(3,indtri1) = inddn-1;
    indtri2 = (m1t+1):2:m2t;
    t(1,indtri2) = inddn;
    t(2,indtri2) = indup+1;
    t(3,indtri2) = indup;    
end

