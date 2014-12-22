function [invalid,mu]=mullwlsk(bw,kernel,xin,yin,win,out1,out2,count)
% use local weighted least squares kernel smoothing.
% input bw(1,2):    specified bandwidth, b(1) for x axis b(2) for y axis
% input kernel:     specified kernel function, candidates:'epan','rectangle'
% input xin(2,n):   matrix of predictors. xin(:,i) is the coordinate of the ith
%                   predictor.
% input yin(n,1):   vector of responses. y(i) is the response valued at x(:,i)
% input win(1,n):   vector of case weights. win(i) is the weight of the ith 
%                   observation.
% input out1(1,m):  the output grid for x-coordinate which is sorted
% input out2(1,k):  the output grid for y-coordinate which is sorted 
% output mu(k,m):   mu(i,j) is the local weighted least square estimate at 
%                   timepoint (out2(i),out1(j))
% deleting observations with weights 0
active=find(win~=0);
xin=xin(:,active);
yin=yin(active);
win=win(active);
invalid=0;
mu=zeros(length(out2), length(out1));
%gap=[];
for i=1:length(out2)
 for j=i:length(out1)
  % locating local window
  if strcmp(kernel,'gauss') ~= 1      %when it is not gaussian kernel, look for the grids that have domain -1 and 1
    list1=find(xin(1,:)>=out1(j)-bw(1)-10^(-6)&xin(1,:)<=out1(j)+bw(1)+10^(-6));
    list2=find(xin(2,:)>=out2(i)-bw(2)-10^(-6)&xin(2,:)<=out2(i)+bw(2)+10^(-6));
    ind=intersect(list1,list2);
  else
     ind = 1:size(xin,2);
  end
  lx=xin(:,ind);
  ly=yin(ind); 
  lw=win(ind);
  % computing weight matrix 
  if length(unique(lx','rows'))>=3
   llx=[(lx(1,:)-out1(j))./bw(1);(lx(2,:)-out2(i))./bw(2)];
   % deciding the kernel used
   k=size(llx,2);
   indd=(1:k)';
   if strcmp(kernel,'epan')
     temp=lw.*(1-llx(1,:).^2).*(1-llx(2,:).^2).*(9/16);
   elseif strcmp(kernel,'rect')
     temp=lw.*ones(1,size(lx,2))/4;
   elseif strcmp(kernel,'gauss')
     temp = lw.*(exp(-0.5*(llx(1,:).^2))/sqrt(2*pi).*(exp(-0.5*(llx(2,:).^2))/sqrt(2*pi)));
   elseif strcmp(kernel,'gausvar')
     temp = lw.*(exp(-0.5*(llx(1,:).^2))/sqrt(2*pi).*(1.25-0.25*llx(1,:).^2).*((exp(-0.5*(llx(2,:).^2))/sqrt(2*pi))).*(1.5-0.5*llx(2,:).^2));
   elseif strcmp(kernel,'quar')
     temp = lw.*((1-llx(1,:).^2).^2).*((1-llx(2,:).^2).^2).*(225/256);
   end   
   
   W=sparse(indd,indd,temp');
   
   % computing design matrix
   X=zeros(length(ly),3);
   X(:,1)=ones(length(ly),1);
   X(:,2)=lx(1,:)'-out1(j);
   X(:,3)=lx(2,:)'-out2(i);
   if nargin == 8
     temp = temp.*count(ind)';
     W1=sparse(indd,indd,temp);
   else
     W1 = W;
   end   
   beta=pinv(X'*W1*X)*X'*W*ly;
   clear X W W1;
   mu(i,j)=beta(1);
   %gap(i,j)=0;
  % if there is only observation in local window, take this value as estimate 
  %elseif length(unique(ly))==1
  % mu(i,j)=unique(ly);
  else
   invalid=1;
   %gap(i,j)=1;
   fprintf(1,'No enough points in local window, please increase bandwidth\n');
   mu = [];
   return;
  end
 end
end
  
if isempty(mu) == 0   
    a = triu(mu,1);    %obtain upper triagular part of the mu matrix
    mu = a+a'+diag(diag(mu)); %assign the lower triangular part of the mu matrix to be the same as the upper triangular part
end

end   
    
    
  
