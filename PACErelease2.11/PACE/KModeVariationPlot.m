%function KModeVariationPlot(yy)
%This function creates the k-th mode of variation plot based on the
%results from FPCA() or FPCder().
%Input yy : returned object from FPCA() or FPCder().
%Input k  : a positive integer (<= no_opt from FPCA(), or FPCder()) 
%           that specifies the k-th mode of variation plot that is 
%           of interest 
%           k can also be a row vector for multiple plots.
%           If missing, then, k = 1:no_opt
%example:
% yy = FPCA(y,t,p);
% KModeVariationPlot(yy)  %k-th mode variation plots for all k from 1 to no_opt
% or 
% p = setDerOptions('nder',0:2);
% yy = FPCder(y,t,p);
% KModeVariationPlot(yy,1) %k-th mode variation plot for k = 1
function KModeVariationPlot(yy, k)
  no_opt = getVal(yy,'no_opt');
  if nargin < 2 
    k = 1:no_opt;
  elseif any(k <= 0)
    fprintf(1,'Error: k must be a positive integer!\n');
    no_opt = [];
    return; 
  elseif any(k > no_opt)
    fprintf(1,['Error: any k cannot be larger than no_opt = ' num2str(no_opt) '\n']);
    no_opt = [];
    return;
  end

  lambda = getVal(yy,'lambda');
  mu = getVal(yy,'mu');
  if iscell(mu)
    mu = mu{1};
  end
  out1 = getVal(yy,'out1');
  phi = getVal(yy,'phi');
  if iscell(phi)
    phi = phi{1};
  end
   for i = 1:length(k)
     figure;
     alpha = [-2 -1 0 1 2].*sqrt(lambda(k(i)));
     cmat = cool(length(alpha));
     for j = 1:length(alpha)
       tmp = mu+alpha(j)*phi(:,k(i))';
       plot(out1, tmp, 'Color',cmat(j,:));
       hold on;  
     end
     xlabel('t');
     ylabel(['\mu(t)+\alpha\phi_' num2str(k(i)) '(t)']);
     title(['K-th mode of variation plot (k = ' num2str(k(i)) ')']);
     legend('\alpha = -2','\alpha = -1','\alpha = 0','\alpha = 1', '\alpha = 2', 'Location','Best')
  end 

end

