function [invalid, xcov, out1, out21, out22]=pc_covNE(y,t,omit,mu,bw_xcov,ngrid,regular,out1,kernel,rcov,xcov)

   tt=cell2mat(t);

   %%%out1=unique(sort(tt));
   a0=min(tt);
   b0=max(tt);
   lint=b0-a0;
   h=lint/(ngrid-1);
   out21=a0:h:b0;
   out22=a0:h:b0;

   y(omit)=[];
   t(omit)=[];

   %rcov = getRawCov(y,t,mu,regular,0);

   %get smoothed covariance surface for x(t) using lwls
   [invalid, xcov]=mullwlsk(bw_xcov,kernel,rcov.tpairn,rcov.cyy',rcov.win,out21,out22);
   if invalid==1
      xcov=[];
   end
