% Example for FPCfam.m 

% Response Y is functional with dense observations;
% Predictor X is functional with dense observations.
% 
% 100 subjects for estimation; 100 subjects for prediction.
% time domain is [0,10] for both X and Y.
% the number of observations for each subject is 
% selected uniformly from {30,...,40}.

% Both functional addtive model and functional linear model 
% are fitted.


p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

rand('twister',sum(10000*clock));
ncohort = 200; npred = 100;
lintx = 10; linty = 10;

%===============Dense data: uniformly samples points======================

lambdax=[4,1]; sigmax = 0.5; sigmay = sqrt(0.1);
pcx = randn(ncohort,2)*diag(sqrt(lambdax));      % 2 FPC score for X
pcy = sum(pcx.^2-repmat(lambdax,ncohort,1),2);

nptx = mysample(30:40,ncohort,1);
npty = mysample(30:40,ncohort-npred,1);
x = cell(1,ncohort); t_x = cell(1,ncohort);
y = cell(1,ncohort-npred); t_y = cell(1,ncohort);
Ey = cell(1,npred); 

% generate the random sparse smaples for estimation
for i = 1:(ncohort-npred)
    t_x{i} = sort(lintx*rand(1,nptx(i)));
    t_y{i} = sort(linty*rand(1,npty(i)));
    %generate the repeated measurements for X with measurement errors N(0,sigmax^2)
    x{i} = mu_true(t_x{i},lintx)+pcx(i,:)*xeig(t_x{i},lintx,length(lambdax))+sigmax*randn(1,nptx(i));
    % calculate the response trajectory Y
    y{i} = mu_true(t_y{i},linty)+pcy(i,:)*xeig(t_y{i},linty,size(pcy,2))+sigmay*randn(1,npty(i)); % zhu: why using xeig? for y? shouldn't it be yeig?
end

% generate the random sparse samples for prediction
ttx = cell2mat(t_x(1:(ncohort-npred)));
tty = cell2mat(t_y(1:(ncohort-npred)));
ntp_y = 40;            % number of timepoints for predicted Y
t_y((ncohort-npred+1):ncohort) = mat2cell(repmat(min(tty):range(tty)/ntp_y:max(tty),npred,1),ones(1,npred))';
for i = (ncohort-npred+1):ncohort
    t_x{i} = min(ttx)+range(ttx)*rand(1,nptx(i));
    x{i} = mu_true(t_x{i},max(ttx))+pcx(i,:)*xeig(t_x{i},max(ttx),length(lambdax))+sigmax*randn(1,nptx(i));
    Ey{i-(ncohort-npred)} = mu_true(t_y{i},max(tty))+pcy(i,:)*xeig(t_y{i},max(tty),size(pcy,2));
end

%===============================================================================

time=cputime;

param_x = setOptions('selection_k','FVE','numBins',50,'regular',0);
param_y = setOptions('selection_k','FVE','numBins',50,'regular',0);
GCV = 1;

% Local polynomial regression of y on pcx_est
[res,xx,yy] = FPCfam(x,t_x,y,t_y,param_x,param_y,[],[],GCV,'gauss',npred);

% Linear regression of y on pcx_est
[res1,xx,yy] = FPCreg(x,t_x,y,t_y,param_x,param_y,0,[],[],npred);

time = cputime-time;
display(['Total time in FPCdiag is : ' num2str(time) ' sec.']);

%=================Plot true and predicted response curves===============================

%plot 9 randomly selected curves
figure;
new_ty = t_y((ncohort-npred+1):ncohort);
newy = getVal(res,'newy');           % estimated response curves
newy1 = getVal(res1,'newy');
id = sort(mysample(1:npred,9,0));
for i = 1:9
    subplot(3,3,i)
    plot(new_ty{id(i)},Ey{id(i)},'k-',new_ty{id(i)},newy{id(i)},'r--',new_ty{id(i)},newy1{id(i)},'b-.');
    title(['Subject ' num2str(id(i))]);
    xlabel('t');
    ylabel('Y(t)');
    legend('true','FAM','LIN','Location','SouthEast')
end

    
    