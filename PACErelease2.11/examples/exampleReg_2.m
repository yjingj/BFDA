% Example for FPCreg.m -- Dense regular case
%
% Random design, 80 subjects, 21 measurements on each subject.
% 50 for estimation and 30 for prediction.
% time domain is [0,10] for predictor X  and [0, 5] for response Y.
% all the observations are taken equidistantly.

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

rand('twister',sum(10000*clock));
ncohort = 80;
npred = 30;
nobs = 20;
lint_x = 10; lint_y = 5;
numeig_x = 3; numeig_y = 3;

tx = 0:lint_x/nobs:lint_x;           % observation time points for X
ty = 0:lint_y/nobs:lint_y;           % observation time points for Y
phi_x = xeig(tx,lint_x,numeig_x);    % 3 eigenfunctions for X
phi_y = xeig(ty,lint_y,numeig_y);    % 3 eigenfunctions for Y          
lambda_x = [3 2 1].^2;
pc_x = randn(ncohort,3)*diag(sqrt(lambda_x));         % 3 FPC score for X
mu_x = mu_true(tx,lint_x);          % true mean function for X
mu_y = mu_true(ty,lint_y);          % true mean function for Y

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

x = cell(1,ncohort); 
t_x = mat2cell(repmat(tx,ncohort,1),ones(1,ncohort))';
y = cell(1,ncohort-npred); 
t_y = mat2cell(repmat(ty,ncohort,1),ones(1,ncohort))';
Ey = cell(1,ncohort);

for i=1:ncohort
    %generate the repeated measurements for X with measurement errors N(0,0.25)
    Ey{i} = mu_y+(pc_x(i,:)*b)*phi_y;
    x{i} = mu_x+pc_x(i,:)*phi_x+0.5*randn(1,length(tx));
    %generate the repeated measurements for Y with measurement errors N(0,0.25)    
    if i <= (ncohort-npred)
        eta = randn(1,numeig_y).*sqrt(lambda_eta);
        y{i} = Ey{i}+eta*phi_y+0.5*randn(1,length(ty));
    end
end

% set parmeter options for FPCA
param_X = setOptions('selection_k','AIC_R','regular',2,'verbose','on');
param_Y = setOptions('selection_k','BIC1','regular',2,'verbose','on');
FIT = 0;
K_x = []; K_y = [];
nsim = 1000;
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
% [res xx yy h d rej pv cb grid_x grid_y] = FPCdiag(x, t_x, y, t_y, param_X, param_Y, FIT, K_x, K_y, npred, bootstrap,nsim,alpha);


time = cputime-time;
display(['Total time in FPCdiag is : ' num2str(time) ' sec.']);

%===========Plot true beta function and estimated beta function===============================
if FIT ==0
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
%plot 9 randomly selected curves
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

