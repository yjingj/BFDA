%Choose the best bandwidth for the mean function through cross-validation method
function bopt=cvfda_lwls(x,t,kernel,nwe,npoly,nder,regular, verbose, bwmuLocal)
if nargin < 9
  bwmuLocal = 0;
end
ncohort=size(t,2);
tt = cell2mat(t);
xx = cell2mat(x);
ind=[];
for i=1:ncohort
    ind=[ind i*ones(1,length(t{i}))];
end
ttn=tt;
xxn=xx;
a0=min(tt);
b0=max(tt);
rang = b0-a0;
dstar = minb(tt, npoly+2);


if regular ~= 2
   h0 = 2.5*dstar;
else
   h0 = dstar;
end

if h0 > rang/4
  h0 = h0*.75;
  display(['Warning: the min bandwidth choice is too big, reduce to ' num2str(h0) '!']);
end

% creat 11 bandwidth candidates
nbw=11;
bw=zeros(1,nbw-1);
%rang=max(tt)-min(tt);
n=length(unique(tt));
for i=1:(nbw-1)
    bw(i)=2.5*rang/n*(n/2.5)^((i-1)/(nbw-1));
end
bw = bw-min(bw)+h0;

ave=zeros(1,length(t{1}));
if regular==2
    for i=1:ncohort
        ave=ave+t{i}/ncohort;
    end
end
    
cv=[];
count=[];
for j=1:(nbw-1);
    time2=cputime;
    cv(j)=0;
    count(j)=0;
    for i=1:ncohort
        out=tt(ind==i);
        obs=xx(ind==i);
        win=ones(1,length(tt));
        win(ind==i)=0;        
        if regular==2 
            xxn=(ave*ncohort-t{i})/(ncohort-1);
            ttn=t{1};
            win=[];
            win=ones(1,length(t{1}));
        end

        [invalid,mu]=lwls(bw(j),kernel,nwe,npoly,nder,ttn,xxn',win,out, bwmuLocal);       

        if invalid==0
           cv(j)=cv(j)+(obs-mu)*(obs-mu)';
           count(j)=count(j)+1;
        end
    end
    time2=cputime-time2;
end
cv=cv(find(count/ncohort>0.90));
bw=bw(find(count/ncohort>0.90));
bopt=bw(find(cv==min(cv), 1, 'first'));
if strcmp(verbose, 'on') == 1
  fprintf(1,['CV bandwidth choice of mean function: ' num2str(bopt) '\n']);           
end

