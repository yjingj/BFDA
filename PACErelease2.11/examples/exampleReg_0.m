% Example for FPCreg.m -- Sparse irregular case
%
% Random design, 300 subjects, 200 for estimation and 100 for prediction. 
% 2~5 measurements on each subject.
% The number of measurements for each subject is 
% equal likely dsitributed among 2,3,4,5.
% time domain is [0,10] for predictor X  and [0,5] for response Y
% The timepoints for each subject are uniformly distributed in the domain.

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

rand('twister',sum(10000*clock));
ncohort = 300;
npred = 100;
lint_x = 10; lint_y = 5;
numeig_x = 3; numeig_y = 3;

ngrid = 50;
tx = 0:lint_x/ngrid:lint_x;          
ty = 0:lint_y/ngrid:lint_y;          
phi_x = xeig(tx,lint_x,numeig_x);      % 3 eigenfunctions for X
phi_y = xeig(ty,lint_y,numeig_y);      % 3 eigenfunctions for Y 
mu_x = mu_true(tx,lint_x);          % true mean function for X
mu_y = mu_true(ty,lint_y);          % true mean function for Y
lambda_x = [3 2 1].^2;
pc_x = randn(ncohort,3)*diag(sqrt(lambda_x));   % 3 FPC score for X

% calculate the true beta function
BETA = zeros(length(tx),length(ty));
b=[1,-0.85,-0.5
   -0.9,0.6,-0.78
   0.79,0.82,-0.39]';
for i=1:numeig_x
    for j=1:numeig_y
        BETA = BETA+b(i,j)*phi_x(i,:)'*phi_y(j,:);
    end
end

lambda_eta = [2 1 0.5];

x = cell(1,ncohort); t_x = cell(1,ncohort);
y = cell(1,ncohort-npred); t_y = cell(1,ncohort);
Ey = cell(1,ncohort);

% calculate the true response E(Y|X)
for i=1:ncohort
    Ey{i} = mu_y+(pc_x(i,:)*b)*phi_y;
end

% generate the random sparse smaples for estimation
for i=1:(ncohort-npred)
    ntp_x = mysample(2:5,1,0);
    t_x{i} = lint_x*rand(1,ntp_x);     % observation timepoints of X_i(not sorted)
    %generate the repeated measurements for X with measurement errors N(0,0.25)
    x{i} = mu_true(t_x{i},lint_x)+pc_x(i,:)*xeig(t_x{i},lint_x,numeig_x)+0.5*randn(1,ntp_x);
    %generate the repeated measurements for Y with measurement errors N(0,0.25)     
    ntp_y = mysample(2:5,1,0);
    t_y{i} = lint_y*rand(1,ntp_y);     % observation timepoints for Y_i(not sorted)
    eta = randn(1,numeig_y).*sqrt(lambda_eta);
    y{i} = interp1(ty,Ey{i},t_y{i},'spline')+eta*xeig(t_y{i},lint_y,numeig_y)+0.5*randn(1,ntp_y);
end
    
% generate the random sparse samples for prediction
ttx = cell2mat(t_x(1:(ncohort-npred)));
tty = cell2mat(t_y(1:(ncohort-npred)));
ntp_y = 50;            % number of timepoints for predicted Y
t_pred = min(tty):range(tty)/ntp_y:max(tty);
for i=(ncohort-npred+1):ncohort
    ntp_x = mysample(2:5,1,0);
    t_x{i} = min(ttx)+range(ttx)*rand(1,ntp_x);
    x{i} = mu_true(t_x{i},max(ttx))+pc_x(i,:)*xeig(t_x{i},max(ttx),numeig_x)+0.5*randn(1,ntp_x);
    t_y{i} = t_pred;
end
    
% set parmeter options for FPCA
param_X = setOptions('selection_k','AIC_R','verbose','on');
param_Y = setOptions('selection_k','FVE','FVE_threshold',0.9,'verbose','on');
FIT = 0;
K_x = []; K_y = [];
nsim = 500;
alpha = 0.05;

time=cputime;
%
%There are two options to run functional regression with diagnostics
%
%Option 1: Run functional regression and its diagnostics seperately
[res, xx, yy] = FPCreg(x, t_x, y, t_y,param_X, param_Y,FIT,K_x,K_y,npred,alpha); %Functional regression
[rej pv] = Rtest(res, xx, yy, nsim, alpha);      % hypothesis testing of Q=0
[h d] =  getDiag(res, xx, yy, param_Y, [], []);  % Residual plots, leverage plot and Cook's distance plot
[cb grid_x grid_y]= getCB(res, xx, yy, [],[]); % 95% Confidence bands plot for Regression function Beta(s,t)

%Option 2: Driver function for FPCreg, getDiag, getCB and Rtest
% bootstrap = 'post'; %or 'pre'
% [res xx yy h d rej pv cb grid_x grid_y] = FPCdiag(x,t_x,y,t_y,param_X,param_Y,FIT,K_x,K_y,npred,bootstrap,nsim,alpha);

time = cputime-time;
display(['Total time in FPCdiag is : ' num2str(time) ' sec.']);

%===========Plot true beta function and estimated beta function===============================
if FIT == 0
    figure;
    subplot(1,2,1)
    mesh(tx,ty,BETA');
    title('True beta function');
    subplot(1,2,2)
    BETAest = getVal(res,'BETA');
    mesh(getVal(BETAest,'grid_x'),getVal(BETAest,'grid_y'),getVal(BETAest,'beta')');
    title('Estimated beta function');
end

%===========Plot true and predicted response curves===================================
% plot 9 randomly selected curves
figure;
newEy = Ey((end-npred+1):end);       % true response curves for prediction part
new_ty = t_y((end-npred+1):end);     % time points corresponding to the newy and newEy
newy = getVal(res,'newy');           % estimated response curves
newcb = getVal(res,'newcb');
id = mysample(1:npred,9,0);
for i = 1:9
    subplot(3,3,i)
    plot(new_ty{id(i)},newEy{id(i)},'k-',new_ty{id(i)},newy{id(i)},'r--',new_ty{id(i)},newcb{1,id(i)},'-.',new_ty{id(i)},newcb{2,id(i)},'-.');
    title(['Subject ' num2str(id(i))]);
    xlabel('t');
    ylabel('Y(t)');
    legend('true','pred','confidence bands', 'Location','SouthEast')
end


