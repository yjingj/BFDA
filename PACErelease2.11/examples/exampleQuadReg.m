% Example for FPCQuadReg.m -- Response Y is continuous scalar
%
% Random design, 80 subjects, 21 measurements on each subject.
% 50 for estimation and 30 for prediction.
% time domain is [0,10] for predictor X and the observations for X are taken equidistantly.

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

rand('twister',sum(10000*clock));
ncohort = 80;
npred = 30;
nobs = 20;
lint_x = 10;

tx = 0:lint_x/nobs:lint_x;          % observation time points for X
numeig_x = 2;
phi_x = xeig(tx,lint_x,numeig_x);           % 2 eigenfunctions for X
pc_x = randn(ncohort,2)*diag([3 2]);      % 2 FPC score for X
mu_x = mu_true(tx,lint_x);          % true mean function for X

x = cell(1,ncohort); 
t_x = mat2cell(repmat(tx,ncohort,1),ones(1,ncohort))';
y = zeros(1,ncohort-npred);

for i=1:ncohort   
    %generate the repeated measurements for X with measurement errors N(0,1)
    true_x = mu_x+pc_x(i,:)*phi_x;
    x{i} = true_x+1*randn(1,length(tx));
    %generate the repeated measurements for Y with measurement errors N(0,1)    
    if i<=(ncohort-npred)
        y(i) = pc_x(i,1)+pc_x(i,2)+pc_x(i,1)^2+pc_x(i,2)^2+pc_x(i,1)*pc_x(i,2)+1*randn(1);
    end
end

% set parmeter options for FPCQuadReg
param_X = setOptions('selection_k','FVE','FVE_threshold',0.9,'error',0,'regular',2,'verbose','on');
K=[];
isNewSub=npred;
[res, xx] = FPCQuadReg(x, t_x, y, param_X, K, isNewSub);
K=getVal(res,'K');
time=cputime;
time = cputime-time;
display(['Total time in FPCdiag is : ' num2str(time) ' sec.']);

%===========Plot true beta and gamma functions and estimated beta and gamma functions===============================
coeffFun = getVal(res,'coeffFun');
phi_est=getVal(xx,'phi');
figure;
BETA=[1,1]*phi_x;
plot(tx,BETA,'k-');
title('Beta function');
hold on;
plot(getVal(coeffFun,'grid'),getVal(coeffFun,'betafun'),'r-');
hold off;
legend('true','estimated');
figure;
GAMMA=phi_x(1,:)'*phi_x(1,:)+phi_x(1,:)'*phi_x(2,:)+phi_x(2,:)'*phi_x(2,:);
mesh(tx,tx,GAMMA);
title('true Gamma function');
figure;
mesh(getVal(coeffFun,'grid'),getVal(coeffFun,'grid'),getVal(coeffFun,'gammafun'));
title('estimated Gamma function');

%===========Plot true and predicted responses==========================================
figure;
newEy = Ey((end-npred+1):end);       % true response curves for prediction part
newy = getVal(res,'newy');           % estimated response curves
plot(newy,newEy,'k.',[min(newy)-0.5,max(newy)+0.5],[min(newy)-0.5,max(newy)+0.5],'r--')
xlabel('Fitted Value'); ylabel('Observed Value')
title('Observed against Fitted Values');
