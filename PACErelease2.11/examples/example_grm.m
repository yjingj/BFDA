% Example for FPCgrm.m
%
% Random design, 200 subjects, 2~5 measurements on each subject.
% Time interval is [0,10].
% The number of measurements for each subject is 
% equal likely dsitributed among 2,3,4,5. The timepoints for each subject 
% are uniformly distributed on [0,10].

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

rand('twister',sum(10000*clock));
ncohort = 200;
lint = 10; ngrid = 50;
tx = 0:lint/ngrid:lint;

noeig = 2; 
lambda = [9 4]; 
phit = xeig(tx,lint,noeig);
muxt = muSimu(tx,lint); 
xi = randn(ncohort,noeig)*diag(sqrt(lambda));

p = setOptions('kernel','gauss','selection_k','BIC1');
muyt = glink(muxt,p.family,p.param);

x = {}; t = {};
ey = {}; y = {};
for i=1:ncohort
    ntp = mysample(2:5,1,0);
    t{i} = lint*rand(1,ntp);    
    phiti = interp1(tx',phit',t{i}')';
    x{i} = interp1(tx,muxt,t{i})+xi(i,:)*phiti;
    ey{i} = glink(x{i},p.family,p.param);
    y{i} = random(p.family,1,ey{i});
end


time = cputime;

[yy] = FPCgrm(y,t,p);

time = cputime-time;
display(['Total time in FPCdiag is : ' num2str(time) ' sec.']);

out1 = getVal(yy,'out1');        %vector of time points for mu, phi and ypred
ymu = getVal(yy,'ymu');          %estimated mean function for observed data Y
xmu = getVal(yy,'xmu');          %estimated mean function for latent Guassian process X
out21 = getVal(yy,'out21');      %vector of time points for covariance
xcovfit = getVal(yy,'xcovfit');  %fitted covariance surface
ycov = getVal(yy,'ycov');        %estimated smooth covariance evaluated at out21 for Y
xcov = getVal(yy,'xcov');        %estimated smooth covariance evaluated at out21 for X
xphi = getVal(yy,'xphi');        %estimated eigenfunctions
no_opt = getVal(yy,'no_opt');    %number of optimal FPC scores
x_pred = getVal(yy,'x_pred');
y_pred = getVal(yy,'y_pred');

figure
subplot(3,2,1)
plot(cell2mat(t),cell2mat(x),'g.');
hold on;
plot(out1,xmu,'r');
xmu_true = muSimu(out1,lint);
plot(out1,xmu_true)
xlabel('t');
ylabel('\mu(t)');
legend('data','fitted','true','Location','Best');
title('Mean function of X');

subplot(3,2,2)
plot(cell2mat(t),cell2mat(y),'y.');
hold on;
plot(out1,ymu,'r');
plot(out1,glink(xmu_true,p.family,p.param))
xlabel('t');
ylabel('\mu(t)');
legend('data','fitted','true','Location','Best');
title('Mean function of Y');

subplot(3,2,3)
xcov_true = zeros(length(out1),length(out1));
phi_true = xeig(out1,lint,noeig);
for i = 1:noeig
    xcov_true = xcov_true+lambda(i)*phi_true(i,:)'*phi_true(i,:);
end
[a b] = meshgrid(out1,out1);
mesh(a,b,xcov_true);
title('True covariance for X');
xlabel('t');
ylabel('t');

subplot(3,2,4)
mesh(out21,out21,xcovfit);
xlabel('t');
ylabel('t');
title('Fitted covariance surface for X');

subplot(3,2,5)
[phi_sign] = checkPhiSign(out1, xphi(:,1), phi_true(1,:)');
plot(out1,phi_sign*xphi(:,1),'r');
hold on
plot(out1,phi_true(1,:), 'b');
xlabel('t')
ylabel('\phi_1(t)')
legend('fitted','true','Location','Best');
title('\phi_1(t) Vs t');

subplot(3,2,6)
plot(out1,phi_true(2,:), 'b');
if no_opt > 1
    hold on
    [phi_sign] = checkPhiSign(out1, xphi(:,2), phi_true(2,:)');
    plot(out1,phi_sign*xphi(:,2),'r');
end
xlabel('t')
ylabel('\phi_2(t)')
title('\phi_2(t) Vs t');
if no_opt > 1
    legend('true','fitted','Location','Best');
else
    legend('true','Location','Best');
end

%plot 9 randomly selected curves
id = mysample(1:ncohort,9,0);   %randomly sample 9 of the curves (without replacement) for plottingfigure

figure
xtrue = cell(1,ncohort);
for i = 1:ncohort
    xtrue{i} = muSimu(out1,lint)+xi(i,:)*phi_true;
end
for k = 1:9
    subplot(3,3,k)
    plot(t{id(k)},x{id(k)},'k*',out1,xtrue{id(k)},'b-',out1,x_pred{id(k)},'r-');
    title(['Subject ' num2str(id(k))]);
    xlabel('t');
    ylabel('X(t)');
    legend('obs','true','pred','Location','Best')
end

figure
for k = 1:9
    subplot(3,3,k)
    plot(t{id(k)},y{id(k)},'k*',out1,y_pred{id(k)},'r-',out1,0.5*ones(1,length(out1)),'b--');
    title(['Subject ' num2str(id(k))]);
    xlabel('t');
    ylabel('X(t)');
    legend('obs','pred','Location','Best')
end

