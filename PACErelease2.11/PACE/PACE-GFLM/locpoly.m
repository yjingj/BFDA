% function for one-dimensional Local Weighted Least Square fit
% input bw                   specified bandwidth. 
% input kernel               kernel function to be used:
%                            'epan' --> epanechnikov: (1-x^2) on [-1,1] 
%                            'rect' --> rectangle:    1 on [-1,1]
%                            'gauss'--> gaussian:     exp(-x^2/2)/sqrt(2*pi) on [-4,4]
%                            'optp' --> opt poylnomial: (1-x^2)^(nwe-1) on [-1,1]
%                            'quar' --> quartic kernel
%                            'gausvar' --> variant of Gaussian kernel
% nwe                        degree of opt polynomial 
% input npoly                degree of polynomial, note order of poly = degree of poly +1
% input nder                 order of derivative  
%                            npoly-nder=odd >0 
% input xin(1,n)             vector of of x-coordinate (predictor)
% input yin(n,1)             vector of y-coordinates of input data.
%                            yin(i) is j-th cohort evaluated at xin(i).
% input win(1,n)             vector of case weights of input data.
% input xou(m)               vector of output grid x-coordinate. Must be Sorted.
% input bwmuLocal:           1    local bandwidth choice is used
%                            0    global bandwidth choice is used (currently always set to 0)
% output invalid             Indicator for error message when bw specified is too small
% output mu(m)               vector of cross mean estimates
% Pseudo-code: (0) Validate setup
%              (1) define the multiple of bandwidth
%              (2) LWLS with different weight functions

function [invalid,mu] = locpoly(bw,kernel,nwe,npoly,nder,xin,yin,win,xou,bwmuLocal)

if nargin < 10
    bwmuLocal = 0;
end
if (npoly < nder) 
    error('Degree of Polynomial should be no less than the order of derivative'); 
    return;
end

actobs=find(win~= 0);
xin=xin(actobs);
yin=yin(actobs);
win=win(actobs);
invalid=0;
% define the multiple of bandwidth
aa=1;
if nwe==-1
    aa=4;
end
mu=zeros(1,length(xou));
gap=[];

if bw > 0
    if bwmuLocal == 1            % local bandwidths
        bw = search_ht(xou,bw);  % this refers to the h0 in the function
    else                         % global bandwidths
        bw = ones(1,length(xou))*bw;
    end
else
    error('Bandwidth choice for mu(t) and/or its derivative must be positive!');
end

% LWLS with different weight functions
for i=1:length(xou)
    % (3-1) Locating local window
    if strcmp(kernel, 'gauss') == 0 && strcmp(kernel, 'gausvar') == 0
        list1=find((xin<=xou(i)+aa*bw(i))); %&(xin~=xou(i)));
        list2=find((xin>=xou(i)-aa*bw(i))); %&(xin~=xou(i)));
        list=intersect(list1,list2);
    else
        list = 1:length(xin);
    end

    lx=xin(list);
    ly=yin(list);
    lw=win(list);

    if length(unique(lx)) >= (npoly+1)
        % Specify weight matrix
        llx=(lx-xou(i))/bw(i);
        if strcmp(kernel,'epan') == 1
            w=lw.*(1-llx.^2)*0.75;
        elseif strcmp(kernel, 'rect') == 1
            w=lw.*ones(1,size(lx,1));
        elseif strcmp(kernel, 'optp') == 1
            w=lw.*(1-llx.^2).^(nwe-1);
        elseif strcmp(kernel,'gauss') == 1
            w=lw.*exp(-0.5*(llx.^2))./sqrt(2*pi);
        elseif strcmp(kernel,'gausvar') == 1
            w=lw.*exp(-0.5*(llx.^2))./sqrt(2*pi).*(1.25-0.25*llx.^2);
        elseif strcmp(kernel, 'quar') == 1
            w = lw.*((1-llx.^2).^2)*(15/16);
        else
            % fprintf(1,'Invalid kernel, Epanechnikov kernel is used!\n');
            w=lw.*(1-llx.^2)*0.75;
        end
        W=diag(w);
        % temp=(1:length(w))';
        % W=sparse(temp,temp,w');
        % Define design matrix
        dx=ones(length(lx),npoly+1);
        for j=1:npoly
            dx(:,j+1)=(xou(i)-lx)'.^j;
        end
        % dxx=sqrt(W)*dx;
        % lyy=sqrt(W)*ly;
        % p=dxx\lyy;
        p=pinv(dx'*W*dx)*dx'*W*ly;
        % Find estimate
        mu(i)=p(nder+1)*factorial(nder)*((-1)^nder);
        gap(i)=0;
        % elseif npoly==1&length(unique(ly))==1&nder==0
        % mu(i)=unique(ly);
    else
        gap(i)=1;
        invalid=1;
        % return;
    end
end
indx=find(gap==0);
if length(indx)>=0.9*length(xou)&length(indx)<length(xou)
    xou1=xou(indx);
    mu1=mu(indx);
    [xou2,id1]=unique(xou1);
    mu2=mu1(id1);
    mu=interp1(xou2,mu2,xou,'spline');
elseif  length(indx)<0.9*length(xou)
    mu='too many gaps, please increase bandwidth';
    invalid = 1;
end

end

