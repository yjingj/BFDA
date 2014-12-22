function [invalid, sigma,  xvar, yvar]=pc_covE(t,out1,bw_xcov,ngrid,cut,kernel,rcov, npoly)

if nargin < 8
  npoly = 1;
end
tt=cell2mat(t);
%%%out1=unique(sort(tt));
a0=min(out1);
b0=max(out1);
lint=b0-a0;
h=range(out1)/(ngrid-1);
out21=linspace(a0,b0,ngrid);
out22=out21;

%rcov = getRawCov(y,t,mu,regular,0);
tpairn = rcov.tpairn;
%get smoothed covariance surface for x(t) using lwls

tneq=find(tpairn(1,:)~=tpairn(2,:));
cyy = rcov.cyy;

if isempty(rcov.count) == 0
 cyy = cyy./rcov.count';    %This is for the case when regular = 1, the raw covariance
end                         %matrix needs to be divided by the number of individual sums
                            %for each element in the matrix. In the case of regular = 2,
                            %the division is n for each of the element. 
cxx=cyy(tneq);
clear rcov;
win1=ones(1,length(cxx));
%[invalid, xcov]=mullwlsk(bw_xcov,kernel,tpairn(:,tneq),cxx',win1,out21,out22);



   %get smoothed variance function for y(t) using lwls
   teq=find(tpairn(1,:)== tpairn(2,:)); 
   vyy=cyy(teq);
   win2=ones(1,length(vyy));
   %yvar is the variance function
   [invalid yvar]=lwls(bw_xcov(1),kernel,1,npoly,0,tpairn(1,teq),vyy',win2,out21,0);

   if invalid==0 
      %estimate variance of measurement error term
      %xvar=diag(xcov)';
      %use quadratic form on diagonal to estimate Var(x(t)) 
      [invalid,xvar]=rotate_mlwls(bw_xcov,kernel,tpairn(:,tneq),cxx',win1,[out21;out22], npoly);

      if invalid==0
         if cut==0
            %sigma=trapzoid(yvar-xvar,out21)/lint; 
            %sigma = romb2(out21,yvar-xvar)/lint;
            sigma = trapz(out21,yvar-xvar)/lint;
         elseif cut==1
            a=a0+lint*0.25;
            b=a0+lint*0.75;
            ind1=find(out21>a & out21<b);
            yvar1=yvar(ind1);
            xvar1=xvar(ind1);
            %sigma=trapzoid(yvar1-xvar1, out21(ind1))*2/lint;
            %sigma = romb(out21(ind1),yvar1-xvar1)*2/lint;
            sigma = trapz(out21(ind1),yvar1-xvar1)*2/lint;
         end
      end
   end
   if sigma < 0
	fprintf(1, 'Warning: estimated sigma is negative, reset to zero now!\n'); 
        sigma = 0;
   end
end
