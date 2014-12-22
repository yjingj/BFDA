% ===========
% Description:
% ===========
% 
%            This program is for aligning warped functional data. It uses
%            a pairwise warping method to obtain the desired alignment (registration)
%            of the random trajectories. The data
%            can be regular (recorded on a regular grid) or sparse. In the
%            latter case, PACE will be used as a preliminary step to
%            evaluate the data on a regular time grid. The algorithm gives
%            the aligned curves and the warping function. Users can also
%            choose to run PACE on the aligned curves.
% 
%            The model is Y(t)=X(h^{-1}(t))+error where Y(t) is the observed noisy
%            curve, X(t) is the aligned curve and h^{-1}(t) is the inverse warping
%            function. The warping function h(t) can be used to align curves, i.e.
%            Y*(t)=Y(h(t)) is output aligned (still with noise, as only the 
%            measurement times are shifted according to the warping).
% 
% 
%            References: Tang, R., M\"{u}ller, H.G. (2008). Pairwise curve
%            synchronization for functional data. Biometrika 95, 875-889.
%            Tang, R., M\"{u}ller, H.G. (2009). Time-synchronized 
%            clustering of gene expression trajectories. Biostatistics 10, 32-45.
%  
% ======
% Usage:
% ======
% 
% function [X,aligned,h,hinv,y_reg,t_reg] = WFPCA(y,t,nknots,lambda,subset,choice,option,p)
% 
% ======
% Input: 
% ======
%      y:          n*m matrix or 1*n cell array where y(i,:) 
%                  or y{i} is a functional trajectory for subject i; PACE
%                  will be used as a preliminary step to evaluate y on a
%                  regular time grid (output t_reg) if y is a cell array.
%      t:          1*m vector (if y is a matrix) or 1*n cell array (if y is
%                  a cell array) of time point grid.
%      nknots:     number of knots used for estimating pairwise warping functions; 
%                  default value is 3.
%      lambda:     penalty parameter used for estimating pairwise warping
%                  functions (larger lambda will force the warping functions
%                  closer to the identity and thus result in less warping);
%                  use 0 as a first default value, which will set the penalty to V*10^-4, 
%                  where V is the average L2 norm of y-mean(y); 
%                  use [] for a second default value (slower option), in which case 
%                  lambda will be selected from the grid of potential values 
%                  V*10^-k (k=1,...,6) as the value which minimizes the variation of y 
%                  after warping (this optimization procedure is time-consuming)
%      subset:     controls the use of subsets to estimate pairwise warping functions.
%                  Then, Pairwise warping functions will only be determined based 
%                  on random subsets of pairs. The smaller the subsets, the faster 
%                  the algorithm  (number of operations is quadratic in subset size), 
%                  while the accuracy of the algorithm declines. 
%                  'Y': use subsets of size min(max(n^0.5,20),n);
%                       the subset is randomly selected for each trajectory [default].
%                  'N': do not use subsets (time-consuming option).
%                  integer: user-specified size of subsets (0<integer<n).
%      choice:     method for estimating warping functions.         
%                  'weighted': use weighted averages of pairwise warping functions and choose 
%                              weights based on inverse pairwise distances  [default].         
%                  'truncated': pairs with large distances (the largest 10%) are truncated and the 
%                               averages of the remaining pairwise distances are used.                  
%      option:     option for running PACE after pairwise warping.
%                  0: perform only pairwise warping [default]. 
%                  1: perform pairwise warping and then PACE.  
%      p:          a structure array obtained from setOptions.m which sets arguments for
%                  PCA.m, not required if option is 0. Default is
%                  setOptions('regular',2,'verbose','off').
% 
% =======
% Output:  
% =======  
%      X:         A cell array containing all returned values from PCA.m; [] if option is 0.  
%      aligned:   n*m matrix where aligned(i,:) is the estimated aligned curve evaluated at
%                 time vector t (for regular data) or t_reg (for sparse data) for the ith subject.  
%      h:         n*m matrix where h(i,:) is the estimated warping function evaluated 
%                 at time vector t (for regular data) or t_reg (for sparse data) 
%                 for the ith subject (see below or Tang and M\"uller, 2008 for the definition).
%      hinv:      n*m matrix where hinv(i,:) is the estimated inverse warping function evaluated at 
%                 time vector t (for regular data) or t_reg (for sparse data) for the ith subject.
%      y_reg:     n*m matrix, if input y is a matrix, the same as y; if y is
%                 a cell array,  the evaluated y on a regular time grid (t_reg)
%                 after PACE is used first and is used for alignment.
%      t_reg:     1*m vector, time grid for y_reg.
% 
%  
% 
%   To get individual values from X, type getVal(X,varname).
%   To obtain the names for X, type names(X).
%   To see an example, check example_WFPCA.m


function [X,aligned,h,hinv,y_reg,t_reg] = WFPCA(y,t,nknots,lambda,subset,choice,option,p)

if nargin<8|isempty(p) p=setOptions('regular',2,'verbose','off'); end
if nargin<7|isempty(option) option=0; end
if nargin<6|isempty(choice) choice='weighted'; end
if nargin<5|isempty(subset) subset='Y'; end         
if nargin<4|isempty(lambda) lambda=0; end            
if nargin<3|isempty(nknots) nknots=3; end

if iscell(y)
    n=length(y);
else
    n=size(y,1);
end
if iscell(t)
    m=30;
else
    m=length(t);
end
% PACE for sparse data
if isa(t,'cell')
    tempt=cell2mat(t);
    t_reg=min(tempt):range(tempt)/(m-1):max(tempt);
    tempp=setOptions('newdata',t_reg,'verbose','on');
    y_input=FPCA(y,t,tempp);  
    tempy=getVal(y_input,'y_pred');
    tempyy=cell2mat(tempy);
    y_reg=vec2mat(tempyy,m);
    clear tempt tempp tempy tempyy;
else
    y_reg=y;
    t_reg=t;
end

% supernorm normalization 
tempsca=max(abs(y_reg'));
if  min(tempsca)~=0
    y_normalized=y_reg.*repmat(1./tempsca',[1 m]); 
else
    for i=1:n
        if tempsca(i)>0
            y_normalized(i,:)=y_reg(i,:)/tempsca(i);
        else
            y_normalized(i,:)=y_reg(i,:);
        end
    end
end
clear tempsca;

for i=1:n
     yfd(i,:)=spapi(2,t_reg,y_normalized(i,:));  % transfer normalized data into spline form
end

for i=1:n
     yfd_old(i,:)=spapi(2,t_reg,y_reg(i,:));   % transfer orignal data into spline form
end

% use subset to estimate warping functions
if strcmp(subset,'Y') 
    subN=floor(min(max(20,sqrt(n)),n)); 
elseif strcmp(subset,'N')
    subN=n;
elseif isa(subset,'numeric')
    subN=subset;
end  

% choose lambda if not specified
if length(lambda)==0 
    Vy=sqrt(sum(trapz(t_reg',(y_normalized-repmat(mean(y_normalized),[n,1]))'.^2))/(n-1));
    Slambda=[];
    for i=1:6
        [tmptX,tmptaligned,tmpth,tmpthinv] = WFPCA(y_reg,t_reg,nknots,Vy*10^-i,'Y',choice,0,[]);
        Slambda(i)=sum(trapz(t_reg',(tmptaligned-repmat(mean(tmptaligned),[n,1]))'.^2));
        clear tmptX tmptaligned tmpth tmpthinv;
    end
    [tmptc,tmptk]=min(Slambda);
    lambda=Vy*10^-tmptk;
    clear tmptc tmptk;
    fprintf(1,['lambda is selected as ' num2str(lambda) '. \n']);
elseif lambda==0
    Vy=sqrt(sum(trapz(t_reg',(y_normalized-repmat(mean(y_normalized),[n,1]))'.^2))/(n-1));
    lambda=Vy*10^-4;
    fprintf(1,['Default lambda value ' num2str(lambda) ' is used. \n']);
end        

indd=1;
hik=[]; 
for k=1:n
     tempind=randperm(n);
     tempind=tempind(1:subN);
     for i=1:subN;
          curvei=yfd(tempind(i),:);
          curvek=yfd(k,:);
          if tempind(i)~=k
              [hik(indd,:),ftemp]=rthik(t_reg,curvei,curvek,nknots,lambda);
              afterdistance(k,i)=mean((rtYH(t_reg,curvei,hik(indd,:))-fnval(curvek,t_reg)).^2);                
              lochat(indd,:)=[k,tempind(i)]; 
              indd=indd+1; 
          else  
              hik(indd,:)=t_reg;
              afterdistance(k,i)=0; 
              lochat(indd,:)=[k,k];  
              indd=indd+1;
          end    
     end
     clear tempind;
end

% take truncated averages or weighted averages of the pairwise warping functions 
h=[];
hinv=[];
aligned=[];
for k=1:n
    weight=[];
    if strcmp(choice,'weighted')
         tempw=afterdistance(k,:); 
         for i=1:subN
             if tempw(i)==0;
                 tempw(i)=inf;
             end
         end
         weight=1./tempw;
         clear tempw;
    elseif strcmp(choice,'truncated')
         weight=afterdistance(k,:)<=quantile(afterdistance(afterdistance>0),.90);            
    end
    if sum(weight)==0
        hinv(k,:)=t_reg;
    else
        weight=weight/sum(weight);
        hinv(k,:)=weight*hik(lochat(:,1)==k,:);
        hinv(k,m)=max(t_reg);
        hinv(k,1)=min(t_reg);
    end
    h(k,:)=fnval(spapi(2,hinv(k,:),t_reg),t_reg);
    h(k,1)=min(t_reg);
    h(k,m)=max(t_reg);
    aligned(k,:)=fnval(yfd_old(k,:),h(k,:));       
end

% tranfer aligned data into cells and run PACE
if option==1
    yy=num2cell(aligned,2)';
    if size(t_reg,1)==m
        tt=repmat(t_reg',[n,1]);
    else
        tt=repmat(t_reg,[n,1]);
    end
    tt=num2cell(tt,2)';
    X=FPCA(yy,tt,p);
else
    X=[];
end

end
    
