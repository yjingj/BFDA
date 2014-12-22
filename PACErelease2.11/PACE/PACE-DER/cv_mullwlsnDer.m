%Choose the best bandwidth choices for the covariance function through Cross-Validation method
function [bw_xcov,crossv]=cv_mullwlsnDer(y,t,mu,ntest,ngrid,regular,error,kernel,npoly, nder, rcov, verbose, bwxcovLocal)
ncohort=size(t,2);      
tt=cell2mat(t);
yy=cell2mat(y);
ind=[];
for i=1:ncohort
    ind=[ind i*ones(1,length(t{i}))];
end
win1=ones(1,length(yy));
out1=unique(sort(tt));
a0=min(out1);
b0=max(out1);
lint=b0-a0;
test=permut(1:ncohort,ntest);
act=0;
for i=1:ntest
   act=act+(length(y{test(i)})>1);
end

%New change:
% Calculate the raw covariance for various cases
% depending on regular,  error
%rcov = getRawCov(y,t,mu,regular, error);
if error == 1
   tpairn = rcov.tpairn;
   tneq=find(tpairn(1,:)~=tpairn(2,:));
   cyy = rcov.cyy;
   rcov.tpairn = tpairn(:,tneq);
   rcov.cxxn=cyy(tneq);
   rcov.win=ones(1,length(rcov.cxxn));
end


%End of New change
dstar = min(range(tt)/4, getMinb(t,out1,regular, npoly));

% create 10 bandwidth candidates
nbw=11;
bw=[];
r=range(tt);
n=length(unique(tt));
for i=1:nbw-1
    bw(i,1)=2.5*r/n*(n/2.5)^((i-1)/(nbw-1));
    bw(i,2)=bw(i,1);
end

bw = bw-min(bw(:,1))+dstar;
if strcmp(kernel,'gauss')
   bw = bw*0.5;
end
%fprintf(1,'bwxcov candidates: \n');
%bw

h=lint/ngrid;
out21=a0:h:b0;
out22=a0:h:b0;
crossv=zeros(1,nbw-1);
count=zeros(1,nbw-1);
len=length(out1);
enter=max(out1(2:len)-out1(1:len-1));
begin=find(bw(:,1)>=enter);
i=begin(1);
continu=1;
while continu==1 
    for j=1:ntest
      omit=test(j);
      if (error==1 & length(y{omit})>1)|error==0
          if regular==0 
              winomit=rcov.win;
              winomit(find(rcov.indx==omit))=0;
              [invalid xcovomit]=mullwlskDer(bw(i,:),kernel,rcov.tpairn,rcov.cxxn',...
                  winomit,npoly, nder, t{omit},t{omit}, bwxcovLocal);
                             
          elseif regular==2 || regular == 1
               winomit=rcov.win;
               if regular == 2
                 [invalid xcovomit]=mullwlskDer(bw(i,:),kernel,rcov.tpairn,rcov.cxxn',...
                  winomit, npoly, nder, t{omit},t{omit}, bwxcovLocal);
               else
                 [invalid xcovomit]=mullwlskDer(bw(i,:),kernel,rcov.tpairn,rcov.cxxn',...
                  winomit, npoly, nder, t{omit},t{omit}, bwxcovLocal,rcov.count);
               end
          end
              
          
          if invalid==0 
            pred=[];
            obs=[];
            k=1;
            for m1=1:length(t{omit})
                for m2=1:length(t{omit})
                    if (error==1 & m1~=m2) | error==0
                        pred=[pred xcovomit(m1,m2)];
                        obs(k)=(y{omit}(m1)-mu(find(out1==t{omit}(m1))))*...
                               (y{omit}(m2)-mu(find(out1==t{omit}(m2))));
                        k=k+1;
                    end
                end
            end
            crossv(i)=crossv(i)+(obs-pred)*(obs-pred)';
            count(i)=count(i)+1;
          end
      end%if
    end%for
    if i==begin(1)
       continu=1;
    elseif error==1 & count(i-1)/act<0.95
       continu=1;
    elseif error==0 & count(i-1)/ntest<0.95
       continu=1;   
    else
       continu=(i<nbw-1)&(crossv(i)<crossv(i-1));
    end
    if continu==1
       i=i+1;
    end
end%while
opt=i-1;
bw_xcov=bw(opt,:);

if strcmp(verbose, 'on') == 1
   fprintf(1,['CV bandwidth choice for COV function: (' num2str(bw_xcov(1)) ',' num2str(bw_xcov(2)) ')\n']);            
end            
        
