function [invalid, xcovDer, phiDer, eigenDer, xcovfitDer]=pc_covNEder(y,t,omit,bw_xcov,lambda, eigen, ngrid,out1,kernel,rcov,bwxcovLocal,npoly, nder, method_int)

   %tt=cell2mat(t);
   out21=linspace(min(out1),max(out1),ngrid);
   out22=out22;

   y(omit)=[];
   t(omit)=[];

   %get smoothed covariance surface for x(t) using lwls
   xcovDer = cell(1,length(nder)-1);
   phiDer = xcovDer;
   eigenDer = xcovDer;
   xcovfitDer = num2cell(zeros(1,length(nder)-1));
   count = rcov.count;
    
   for dg = 2:length(nder)
      if isempty(count) 
       [invalid, xcovDer{dg-1}]=mullwlskDer(bw_xcov(dg,:),kernel{dg},rcov.tpairn,rcov.cyy',rcov.win, npoly(dg), nder(dg), out21,out22,bwxcovLocal(dg));
      else
        [invalid, xcovDer{dg-1}]=mullwlskDer(bw_xcov(dg,:),kernel{dg},rcov.tpairn,rcov.cyy',rcov.win, npoly(dg), nder(dg), out21,out22,bwxcovLocal(dg),count);
      end
   end

   if invalid == 1
       xcovDer = []; phiDer = []; eigenDer = []; xcovfitDer = [];
       return;
   else
       for dg = 2:length(nder)
          for k = 1:size(eigen,2)
             tmp = xcovDer{dg-1}.*repmat(eigen(:,k),1,ngrid);
             if method_int == 1
                tmp = trapz(out21,tmp)/lambda(k);
             elseif method_int == 2
                tmp = romb2(out21,tmp)/lambda(k);
             else
                tmp = gaussquad(out21,tmp)/lambda(k);
             end
             eigenDer{dg-1}(:,k) = tmp;
             phiDer{dg-1}(:,k) = interp1(out21,tmp,out1,'spline');
             xcovfitDer{dg-1} = xcovfitDer{dg-1} + lambda(k)*eigen(:,k)*tmp';
          end
       end
   end
end

