% Example for FSVD.m
%
% Sparse case:
% Random design, 200 subjects. 
% 3~5 measurements on each subject.
% The number of measurements for each subject is 
% equal likely dsitributed among 3,4,5.
% time interval is [0,10] for X  and [0,5] for Y
% The timepoints for each subject are uniformly distributed in the domain.

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

rand('twister',sum(10000*clock));
ncohort = 100;
lint_x = 10; lint_y = 5;
x = cell(1,ncohort); t_x = cell(1,ncohort);
y = cell(1,ncohort); t_y = cell(1,ncohort);

ngrid = 40;
tx = 0:lint_x/ngrid:lint_x;          
ty = 0:lint_y/ngrid:lint_y;          
phi_x = xeig3(tx,lint_x);           
phi_y = xeig3(ty,lint_y);           
phi_x = phi_x([2 3 1],:);
% calculate the true cross-covariance surface
Sig= [8    3   -2    3         0         0;
    3    4    1         0    1.5         0;
   -2    1    3         0         0    0.5000;
    3         0         0    6   -2    1;
         0    1.5         0   -2    4.5   -1.5;
         0         0    0.5    1   -1.5    3.25];
Sig2=Sig^(1/2);
ccov_true = phi_x' * Sig(1:3,4:6) * phi_y;

% generate the random sparse samples for estimation
% pc_x = randn(ncohort,3)*diag([4 2 1]);   % 3 FPC score for X
pc = randn(ncohort,6)*Sig2;
simx=pc(:,1:3)*phi_x;
simy=pc(:,4:6)*phi_y;
for i=1:ncohort
    ntp_x = mysample(3:5,1,0);
    t_x{i} = lint_x*rand(1,ntp_x);     % observation timepoints of X_i(not sorted)
    %generate the repeated measurements for X with measurement errors N(0,0.25)
    x{i} = mu_true(t_x{i},lint_x)+interp1(tx,simx(i,:),t_x{i},'spline')+0.5*randn(1,ntp_x);
    %generate the repeated measurements for Y with measurement errors N(0,0.25)     
    ntp_y = mysample(3:5,1,0);
    t_y{i} = lint_y*rand(1,ntp_y);     % observation timepoints for Y_i(not sorted)
    y{i} = mu_true(t_y{i},lint_y)+interp1(ty,simy(i,:),t_y{i},'spline')+0.5*randn(1,ntp_y);
end
    
time=cputime;
%
%Run functional singular component analysis
res = FSVD(x, t_x, [], y, t_y, [], []);

time= cputime-time;
display(['Total time in FSVD is : ' num2str(time) ' sec.']);


%===========Plot true cross-covariance and estimated cross-covariance===============================
% figure;
% subplot(1,2,1)
% mesh(tx,ty,ccov_true);
% title('True cross-covariance surface');
% ax = axis;
% subplot(1,2,2)
% ccovfit = getVal(res,'ccovfit');
% mesh(getVal(res,'out2x'),getVal(res,'out2y'),ccovfit');
% title('Estimated cross-covariance surface');
% axis(ax);

figure;
subplot(1,3,1)
mesh(tx,ty,ccov_true);
title('True cross-covariance surface');
ax = axis;
subplot(1,3,2)
ccovfit = getVal(res,'ccovfit');
mesh(getVal(res,'out2x'),getVal(res,'out2y'),ccovfit');
title('Estimated cross-covariance surface');
axis(ax);
subplot(1,3,3)
mesh(getVal(res,'out2x'),getVal(res,'out2y'),getVal(res,'ccov_s')');
title('Smoothed cross-covariance surface');
axis(ax);

lambda = getVal(res,'lambda');
fprintf(1, 'The first singular value: %f\n', lambda(1));
% getVal(res,'bwccov')


