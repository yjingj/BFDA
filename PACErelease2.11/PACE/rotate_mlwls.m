function [invalid,mu]=rotate_mlwls(bw,kernel,xin,yin,win,d, npoly)
% input bw:         1*2 vector, bandwidth of smoothing
% input kernel:     specified kernel function, candidates:'epan','rectangle' or 'gauss'
% input xin(2,n):   matrix of predictors. xin(:,i) is the coordinate of the ith
%                   predictor.
% input yin(n,1):   vector of responses. y(i) is the response valued at x(:,i)
% input win(1,n):   vector of case weights. win(i) is the weight of the ith 
%                   observation.
% input d:          2*m, coordinates of output grid, diagonal in this case
% input npoly :   a integer denotes the degree of polynomial to fit
% output mu(k,m):   mu(i,j) is the local weighted least square estimate at 
%                   timepoint (out2(i),out1(j))
% deleting observations with weights 0
if nargin < 7
  npoly = 1;
end

active=find(win~=0);
xin=xin(:,active);
yin=yin(active);
win=win(active);
% rotating coordinates of predictors by pi/4
R=sqrt(2)/2*[1,-1;1,1];
xn=R*xin;
yn=yin;
dn=R*d;
mu=[];
invalid=0;
% minimizing local weigthed least squares
m=size(d,2);
for i=1:m
    % locating local window
    if strcmp(kernel, 'gauss') == 0 && strcmp(kernel,'gausvar') == 0
      list1=find(xn(1,:)>=dn(1,i)-bw(1)&xn(1,:)<=dn(1,i)+bw(1));
      list2=find(xn(2,:)>=dn(2,i)-bw(2)&xn(2,:)<=dn(2,i)+bw(2));
      ind=intersect(list1,list2);
    else
      ind = 1:size(xn,2);
    end
    lx=xn(:,ind);
    ly=yn(ind);
    lw=win(ind);
    % computing weight matrix 
    if length(ly)>=npoly+1
        llx=[(lx(1,:)-dn(1,i))./bw(1);(lx(2,:)-dn(2,i))./bw(2)];
        % deciding the kernel used
       if strcmp(kernel,'epan')
            w=lw.*(1-llx(1,:).^2).*(1-llx(2,:).^2).*(9/16);
 	elseif strcmp(kernel,'rect')
            w=lw.*ones(1,size(lx,2))/4;
 	elseif strcmp(kernel,'gauss')
            w = lw.*(exp(-0.5*(llx(1,:).^2))/sqrt(2*pi).*(exp(-0.5*(llx(2,:).^2))/sqrt(2*pi)));
        elseif strcmp(kernel,'gausvar')
            w = lw.*(exp(-0.5*(llx(1,:).^2))/sqrt(2*pi).*(1.25-0.25*llx(1,:).^2).*((exp(-0.5*(llx(2,:).^2))/sqrt(2*pi))).*(1.5-0.5*llx(2,:).^2));
        elseif strcmp(kernel,'quar')
            w = lw.*((1-llx(1,:).^2).^2).*((1-llx(2,:).^2).^2).*(225/256);
        end
        temp=(1:length(w))';
        W=sparse(temp,temp,w');
        % computing design matrix
        X=[];
        X(:,1)=ones(length(ly),1);
        X(:,2)=(lx(1,:)'-dn(1,i)).^2;
        X(:,3)=lx(2,:)'-dn(2,i);
        %XX=sqrt(W)*X;
        %YY=sqrt(W)*ly;
        %beta=XX\YY;
        beta=pinv(X'*W*X)*X'*W*ly;
        mu(i)=beta(1);
        % if there is only observation in local window, take this value as estimate 
    elseif length(ly)==1
        mu(i)=ly;
    else
        disp('no points in local window, please increase bandwidth');
        invalid=1;
        return;
    end
end
  
