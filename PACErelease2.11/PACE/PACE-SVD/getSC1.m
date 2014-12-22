function [res invalid]=getSC1(x,t_x,xx,y,t_y,yy,bwccov,nsvd,regular)

ngrid=51;
res=[];
sigma1_x = getVal(xx,'rho_opt');
sigma1_y = getVal(yy,'rho_opt');
%obtain smoothed cross covariance
x_mu=getVal(xx,'mucopy');
y_mu=getVal(yy,'mucopy');
out2x=getVal(xx,'out21copy');
out2y=getVal(yy,'out21copy');
[rccov, tpair] = getRawCCOV(x,t_x,y,t_y,x_mu,y_mu,1,regular);
win=ones(1,length(tpair));
if isempty(bwccov)
    bwccov = gcv2_mullwlsn(t_x,t_y,ngrid,regular,'epan',tpair,rccov,win);
    bwccov = bwfactor(t_x, t_y)*bwccov;
end

[invalid, ccov_s]=mullwlsk(bwccov,'epan',tpair,rccov',win,out2x,out2y);
if invalid == 1
    return;
end

ccov=ccov_s;

%Singular components
[sc_x d sc_y] = svds(ccov, ngrid-2);
hx = trapz(out2x, sc_x.^2)';
hy = trapz(out2y,  sc_y.^2)';
sc_x = sc_x / diag(sqrt(hx));
sc_y = sc_y / diag(sqrt(hy));
d = diag(d).*sqrt(hx.*hy);
d = d';
if isempty(nsvd)
   FVE = cumsum(d)./sum(d); 
   nsvd = find(FVE > .75, 1, 'first');
   createSSP(FVE, nsvd, 'X and Y');   
end
lambda = d(1:nsvd);
lambda_x = diag(sc_x' * getVal(xx,'xcovfitcopy') * sc_x);
lambda_y = diag(sc_y' * getVal(yy,'xcovfitcopy') * sc_y);
lambda_x = lambda_x .* hx.^2; 
lambda_y = lambda_y .* hy.^2; 
lambda_x = lambda_x(1:nsvd);
lambda_y = lambda_y(1:nsvd);

out_x=getVal(xx,'out1copy');
out_y=getVal(yy,'out1copy');
ccovfit = sc_x(:,1:nsvd) * diag(lambda) * sc_y(:,1:nsvd)';
sc_x=interp1(out2x, sc_x, out_x, 'spline');
sc_y=interp1(out2y, sc_y, out_y, 'spline');
sc_x = sc_x(:,1:nsvd);
sc_y = sc_y(:,1:nsvd);
% save temp;
for i=1:size(sc_x,2)
    sc_x(:,i) = sc_x(:,i)/sqrt(trapz(out_x,sc_x(:,i).^2));
end
for i=1:size(sc_y,2)
    sc_y(:,i) = sc_y(:,i)/sqrt(trapz(out_y,sc_y(:,i).^2));
end

resnames = {'nsvd' 'sigma1_x' 'sigma1_y' 'xx' 'yy' 'sc_x', 'sc_y', 'out_x', 'out_y', 'out2x', 'out2y' 'ccov_s' 'ccovfit', 'lambda', 'lambda_x', 'lambda_y', 'bwccov', 'rccov', 'tpair','names'};
res = {nsvd, sigma1_x, sigma1_y, xx, yy, sc_x, sc_y, out_x, out_y, out2x, out2y, ccov_s, ccovfit, lambda, lambda_x, lambda_y, bwccov, rccov, tpair, resnames};
end

function createSSP(FVE, no_opt, yname)
if nargin < 3
  yname = 'X and Y';
end

figure;
plot(0:length(FVE), [0 FVE]*100, '--ro', 'LineWidth',2,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
xlabel('\bf{No. of Singluar Components}');
ylabel('\bf{FVE (%)}');
title(['\bf{Fraction of variance explained by No. of SC for function ' yname '}'])
hold on;
plot(linspace(0,no_opt,30), ones(30,1)*FVE(no_opt)*100,'b',...
     ones(30,1)*no_opt, linspace(0,FVE(no_opt)*100,30),'b');
text(no_opt+0.2, FVE(no_opt)*100-10, {[' k = ' num2str(no_opt) ', FVE = ' num2str(roundoff(FVE(no_opt)*100,3)) '%'] ' (final choice)'});
axis([0 length(FVE)+1 0 101]);
hold off;

end

function factor=bwfactor(t1,t2)
    ncohort=length(t1);     % obtain the number of curves or subjects
    ni = zeros(ncohort,2);
    for i = 1:ncohort
        ni(i,1)= length(t1{i});
        ni(i,2)= length(t2{i});
    end
    
    nmean = mean(mean(ni));
    if nmean > 20 %dense
        factor = 1.3;
    else
        factor = interp1([2 3 4 5 10 20], [1 0.9 0.9 0.8 1.1 1.3], nmean);
    end

end