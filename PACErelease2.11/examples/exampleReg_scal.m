% Example for FPCreg.m -- Response Y is continuous scalar
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
numeig_x = 3;
phi_x = xeig(tx,lint_x,numeig_x);           % 3 eigenfunctions for X
pc_x = randn(ncohort,3)*diag([3 2 1]);      % 3 FPC score for X
mu_x = mu_true(tx,lint_x);          % true mean function for X

b=[2 2 1];
BETA = b*phi_x;

x = cell(1,ncohort); 
t_x = mat2cell(repmat(tx,ncohort,1),ones(1,ncohort))';
y = zeros(1,ncohort-npred); 
t_y = [];
Ey = zeros(1,ncohort);

for i=1:ncohort   
    %generate the repeated measurements for X with measurement errors N(0,1)
    true_x = mu_x+pc_x(i,:)*phi_x;
    x{i} = true_x+1*randn(1,length(tx));
    betax = diag(true_x)*BETA';
    Ey(i) = trapz(tx,betax');
    %generate the repeated measurements for Y with measurement errors N(0,1)    
    if i<=(ncohort-npred)
        y(i) = Ey(i)+1*randn(1);
    end
end

% set parmeter options for FPCA
param_X = setOptions('selection_k','AIC_R','regular',2,'verbose','on');
param_Y = [];
FIT = 0;
K_x = []; K_y = [];
nsim = 1000;
alpha = [];

time=cputime;
%
%There are two options to run functional regression with diagnostics
%
%Option 1: Run functional regression and its diagnostics seperately  
[res, xx, yy] = FPCreg(x, t_x, y, t_y,param_X, param_Y,FIT,K_x,K_y,npred,alpha); %Functional regression
[rej pv] = Rtest(res, xx, y, nsim, alpha);      % hypothesis testing of Q=0
[h d] =  getDiag(res, xx, y, param_Y, [], []);  % Residual plots, leverage plot and Cook's distance plot
[cb grid_x grid_y]= getCB(res, xx, y, [],[]); % 95% Confidence bands plot for Regression function Beta(s,t)

%Option 2: Driver function for FPCreg, getDiag, getCB and Rtest
% bootstrap = 'post'; %or 'pre'
% [res xx yy h d rej pv cb grid_x grid_y] = FPCdiag(x,t_x,y,t_y,param_X,param_Y,FIT,K_x,K_y,npred,bootstrap,nsim,alpha);

time = cputime-time;
display(['Total time in FPCdiag is : ' num2str(time) ' sec.']);

%===========Plot true beta function and estimated beta function===============================
figure;
plot(tx,BETA,'k-');
title('Beta function');
hold on;
BETAest = getVal(res,'BETA');
plot(getVal(BETAest,'grid_x'),getVal(BETAest,'beta'),'r-');
hold off;

%===========Plot true and predicted responses==========================================
figure;
newEy = Ey((end-npred+1):end);       % true response curves for prediction part
newy = getVal(res,'newy');           % estimated response curves
plot(newy,newEy,'k.',[min(newy)-0.5,max(newy)+0.5],[min(newy)-0.5,max(newy)+0.5],'r--')
xlabel('Fitted Value'); ylabel('Observed Value')
title('Observed against Fitted Values');
