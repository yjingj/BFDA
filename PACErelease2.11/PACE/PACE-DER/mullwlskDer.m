function [invalid,mu]=mullwlskDer(bw,kernel,xin,yin,win,npoly, nder, out1,out2, bwxcovLocal,counts)
% use local weighted least squares kernel smoothing.
% input bw(1,2):    specified bandwidth, b(1) for x axis b(2) for y axis
% input kernel:     specified kernel function, candidates:'epan','rect' or 'gauss'
% input xin(2,n):   matrix of predictors. xin(:,i) is the coordinate of the ith
%                   predictor.
% input yin(n,1):   vector of responses. y(i) is the response valued at x(:,i)
% input win(1,n):   vector of case weights. win(i) is the weight of the ith 
%                   observation.
% input out1(1,m):  the output grid for x-coordinate which is sorted (derivative)
% input out2(1,k):  the output grid for y-coordinate which is sorted (smoothing)
% bwxcovLocal    :  1 local bandwidth choice 
%                   0 global bandwidth choice   (currently always set to 0)
% counts         :  integer of counts, only used if regular = 1
% output mu(k,m):   mu(i,j) is the local weighted least square estimate at 
%                   timepoint (out2(i),out1(j))
%                   out1 is the derivative direction and out2 is the smoothing direction
%                   Smoothing directions is row-wise. Derivative directions is column-wise.
if npoly < nder 
     error('Degree of Polynomial should be no less than the order of derivative!');
     return;
end

% deleting observations with weights 0
active=find(win~=0);
xin=xin(:,active);
yin=yin(active);
win=win(active);

invalid=0;
mu=zeros(length(out2),length(out1));

%Add by Bitao to allow for local bandwidth choice
if bw > 0
   if bwxcovLocal == 1
     bw1 = search_ht(out1, bw(1));   %this refers to the h0 in the function
     bw2 = search_ht(out2, bw(2));
     bw = [bw1; bw2];
   else
         bw = repmat(bw',1,length(out1));
   end
else
   error('Bandwidth choice for mu(t) must be positive!');
end

jstart = 1; %to speed up computation time if nder = 0
for i=1:length(out2)   %for each row of grids in mu(out2(i),out1(j))
  if nder == 0
    jstart = i;
  end
  for j=jstart:length(out1)  %for each column of row out2(i)
  % locating local window
  if strcmp(kernel, 'gauss') == 0 && strcmp(kernel,'gausvar') == 0
    list1=find(xin(1,:)>=out1(j)-bw(1,j)-10^(-6)&xin(1,:)<=out1(j)+bw(1,j)+10^(-6));
    list2=find(xin(2,:)>=out2(i)-bw(2,i)-10^(-6)&xin(2,:)<=out2(i)+bw(2,i)+10^(-6));
    ind=intersect(list1,list2);
  else
    ind = 1:size(xin,2);
  end
  lx=xin(:,ind);
  ly=yin(ind);
  lw=win(ind);
  % computing weight matrix 
  if length(unique(lx','rows'))>=npoly+2
   llx=[(lx(1,:)-out1(j))./bw(1,j);(lx(2,:)-out2(i))./bw(2,i)];
   % deciding the kernel used
   k=size(llx,2);
   indd=(1:k)';
   if strcmp(kernel,'epan')
     temp=lw.*(1-llx(1,:).^2).*(1-llx(2,:).^2).*(9/16);
     %W=diag(temp);
   elseif strcmp(kernel,'rect')
     temp=lw.*ones(1,size(lx,2))/4;
     %W=diag(temp);
   elseif strcmp(kernel, 'gauss')
     temp = lw.*(exp(-0.5*(llx(1,:).^2))/sqrt(2*pi).*(exp(-0.5*(llx(2,:).^2))/sqrt(2*pi)));
   elseif strcmp(kernel,'gausvar')
     temp = lw.*(exp(-0.5*(llx(1,:).^2))/sqrt(2*pi).*(1.25-0.25*llx(1,:).^2).*((exp(-0.5*(llx(2,:).^2))/sqrt(2*pi))).*(1.5-0.5*llx(2,:).^2));
   elseif strcmp(kernel,'quar')                                                                                                                            
      temp = lw.*((1-llx(1,:).^2).^2).*((1-llx(2,:).^2).^2).*(225/256);                                                                                   
   else
     fprintf(1,['Warning: kernel ' kernel ' is not found! Reset to epan now!\n']);
     temp=lw.*(1-llx(1,:).^2).*(1-llx(2,:).^2).*(9/16);
   end
   
   W = sparse(indd,indd, temp);
   
   if nargin == 11
     temp = temp.*counts(ind)';
     W1=sparse(indd,indd,temp);
   else
     W1 = W;
   end

   % computing design matrix
   X=zeros(length(ly),npoly+2);
   X(:,1)=ones(length(ly),1);
   for jj = 1:npoly
     tmp = lx(1,:)'-out1(j);
     X(:,jj+1)=tmp.^jj;
   end   
   X(:,jj+2)=lx(2,:)'-out2(i);
   %XX=sqrt(W)*X;
   %YY=sqrt(W)*ly;
   %beta=XX\YY;
   beta=pinv(X'*W1*X)*X'*W*ly;
   %mu(i,j)=beta(1);
   mu(i,j) = beta(nder+1)*factorial(nder);
   %gap(i,j)=0;
  % if there is only observation in local window, take this value as estimate 
  %elseif length(unique(ly))==1
  % mu(i,j)=unique(ly);
  else
   %disp('no points in local window, please increase bandwidth');
   invalid=1;
   %gap(i,j)=1;
   mu='no enough points in local window, please increase bandwidth'
   return;
  end
 end
end

if nder == 0 && isempty(mu) == 0
     a = triu(mu,1);           %obtain upper triagular part of the mu matrix
     mu = a+a'+diag(diag(mu)); %assign the lower triangular part of the mu matrix to be the same as the upper triangular part
end

end

function [h] = search_h(t,c1,c2, c3)
    h = c1+c2*exp(-((t-75).^2)./c3);
end
  
   
function [h] = search_ht(t, h0)

   h = h0*(1+exp(-((t-75).^2)./6000));

end   
   
   
    
    
  

