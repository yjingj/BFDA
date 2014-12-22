%function [bw_xcov,gcv]=gcv_mullwlsnDer(t, ngrid,regular,error,kernel,npoly, nder, rcov, verbose, bwxcovLocal)
%This function computes the bandwidth choices for the smooth covariance surface
%using GCV method by pooling all time point pairs from all subjects together.
function [bw_xcov,gcv]=gcv_mullwlsnDer(t, ngrid,regular,error,kernel,npoly, nder, rcov, verbose, bwxcovLocal)
tt = cell2mat(t);
out1=unique(tt);
clear tt;
a0=min(out1);
b0=max(out1);

h0 = getMinb(t,out1,regular, npoly);
if strcmp(kernel,'gauss') || strcmp(kernel,'gausvar')
     %h0 = h0*interp1(1:10,[0.15 0.18 0.2 0.25 0.28 0.33 0.38 0.43 0.46 0.49],npoly,'spline');
     if isempty(h0)
       h0 = max(b0);
     end
     h0 = h0*interp1(1:10,[0.2 0.23 0.28 0.31 0.34 0.37 0.4 0.43 0.46 0.49],npoly,'spline');
end

if isempty(h0)
  fprintf(1,'Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n');                                 
  bw_xcov = [];gcv = [];                                                                                                                                   
  return;       
end

%t1 = cputime;
% Calculate the raw covariance for various cases
% depending on regular, error
%rcov = getRawCov(y,t,mu,regular, error);
clear t;

rcovcount = rcov.count;

if error == 1
  tpairn = rcov.tpairn;
  tneq=find(tpairn(1,:)~=tpairn(2,:));
  cyy = rcov.cyy;
  tpairn = tpairn(:,tneq);
  cxxn=cyy(tneq);
  win=ones(1,length(cxxn));
  if regular == 1
     rcovcount = rcovcount(tneq);
  end
else
  tpairn = rcov.tpairn;
  cxxn = rcov.cxxn;
  win =  rcov.win;
end

clear rcov cyy tneq;
%fprintf(1,['Time after the raw covariance calculation: ' num2str(cputime-t1) '\n']);
% create 10 bandwidth candidates

%Alternative choice: h0= 2.5*dstar; or h0 = 1.1*dstar;
%h0 = dstar;
N = length(cxxn);   %No. of pairs for the covariance
r = range(out1);
clear out1;

%Alternative way to set h0
%dstar = minb(tpairn(1,:),10);
%h0 = sqrt(2.5*dstar/range(tpairn(1,:)));

h0 = min(h0, r);
q = (r/(4*h0))^(1/9);
bw = q.^(0:9).*h0;        %create 10 h candidates
%fprintf(1, 'bwxcov candidates:\n')
bw = sort(repmat(bw',1,2));     %same bandwidth choice for each direction of the covariance

%bw

k0 = mykernel(0, kernel);
if strcmp(kernel, 'gausvar') == 1
  k1 = mykernel(0,'gausvar1');
else
  k1 = k0;
end

out21 = linspace(a0,b0,ngrid); %then length(out21) = ngrid
%bw = sort(bw);
%curTIME = cputime;

leave = 0;
nleave = 0;
tooSparse = 0;
while leave == 0
   gcv = Inf*ones(size(bw,1),1);
   for k = 1:size(bw,1)

  %    fprintf(1,['k = ' num2str(k) '\n']);     
       if regular == 1
          [invalid, xcov]= mullwlskDer(bw(k,:), kernel, tpairn, cxxn', win, npoly, nder, out21, out21, bwxcovLocal, rcovcount);    
       else
          [invalid, xcov]= mullwlskDer(bw(k,:), kernel, tpairn, cxxn', win, npoly, nder, out21, out21, bwxcovLocal);    
       end
%      fprintf(1,'Done with smooth the raw covariance\n');
       %interpolate the smooth covariance from (out21,out21) to (tpairn(1,:), tpairn(2,:))
       if invalid ~= 1
          newxcov = interp2(out21,out21,xcov, tpairn(1,:),tpairn(2,:),'spline');
          clear xcov;         
   %      fprintf(1,'Done with the interpolation\n');
          if regular == 1
             cvsum = (cxxn./(rcovcount')-newxcov)*(cxxn./(rcovcount')-newxcov)';
	  else
             cvsum = (cxxn-newxcov)*(cxxn-newxcov)';
          end
          clear newxcov;      
          bottom = 1-(1/N)*k0*k1*(r/bw(k))^2;
          gcv(k) = cvsum/(bottom)^2;
          tmp = gcv(~isinf(gcv));
          if length(tmp) > 1 && gcv(k) > gcv(k-1)
	     leave = 1;
             break;
          end
       end
   end
   if all(gcv == Inf)
        if nleave == 0 && bw(10,1) < r
             bw_xcov = bw(10,:);
             tooSparse = 1;
        else
             fprintf(1,'Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n');
             bw_xcov = []; gcv = [];
             return;
        end
   else
        bw_xcov = bw(find(gcv == min(gcv),1,'first'),:);
   end

   if bw_xcov(1) == r
        leave = 1;
        fprintf(1,'Warning: data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n');
   elseif bw_xcov(1) == bw(10,1) && nleave == 0
        if (tooSparse == 1) || (sum(gcv == Inf) == 9)
            fprintf(1,'Warning: data is too sparse, retry with larger bandwidths\n');
            h0 = bw(10,1)*1.01;
        else
            fprintf(1, 'Bandwidth candidates are too small, retry with larger choices now!\n');
            h0 = bw(9,1);
        end
        newr = (0.5:0.05:1).*r;
        id = find(h0 < newr, 1, 'first');
        q = (newr(id)/h0)^(1/9);
        bw = q.^(0:9).*h0;
        
        bw = sort(repmat(bw',1,2));     %same bandwidth choice for each direction of the covariance
        if strcmp(verbose, 'on') == 1
          fprintf(1, 'New bwxcov candidates:\n');
          bw
        end
   elseif bw_xcov(1) < bw(10,1) || nleave > 0
        leave = 1;
   end
   nleave = nleave+1;

end  %end of while loop

   if strcmp(kernel,'gauss') == 0 && strcmp(verbose, 'on') == 1
      fprintf(1,['GCV bandwidth choice for COV function (npoly = ' num2str(npoly) '): (' num2str(bw_xcov(1)) ',' num2str(bw_xcov(2)) ')\n']);
   end
end
