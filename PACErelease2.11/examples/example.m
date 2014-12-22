
%  Example file for how to use FPCA.m.
%  Random design, 200 subjects with 1 to 4 measurements on each subject.
%  Time interval is [0,10]. The number of measurements for each subject is
%  uniformly distributed on {1,2,3,4}. The timepoints for each subject
%  are uniformly distributed on [0,10]. 
%
%  In this example, the goal is to predict the trajectories.

%generate data set
%clear all;

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
rand('twister',sum(10000*clock));
mtp = 4;          %at most 4 repeated measurements in the simulated data for the sparse case
%mtp = 20;        %20 repeated measurements in the simulated data for
                  %regular and regular with missing cases
ncohort=200;    %200 subjects in the simulated data
lint=10;
y=cell(1,ncohort);
t=cell(1,ncohort);
newy = y;
newt = t;
xi=[];
regular = 0;

% Case iii) regular data with missing values (regular = 1)   
% ni_missing = zeros(1,ncohort);
% for i = 1:ncohort
%  ni_missing(i) = poissrnd(mysample(1:floor(mtp/3),1,0),1,1); 
%  if ni_missing(i) >= mtp
%     ni_missing(i) = mtp-1;
% end
% end

for i=1:ncohort

%   Case i) Sparse and irregular case (regular = 0) 
   ntp=ceil(mtp*rand(1));
   t{i}=lint*rand(1,ntp);                 
   newt{i} = lint*rand(1,ntp);

   %Case ii) complete balance case (regular = 2) 
%    regular = 2;
%    t{i} = linspace(0,lint,mtp);                      
%    newt{i} = lint*rand(1,ntp);

   %Case iii) regular data with missing values (regular = 1)
%    regular = 1;
%    t{i} = linspace(0,lint,mtp);                    
%    newt{i} = linspace(0,lint,mtp);
%    if ni_missing(i) > 0
%     id_missing = mysample(1:mtp, ni_missing(i),0);
%     t{i}(id_missing) = [];
%     newt{i}(mysample(1:mtp,ni_missing(i),0)) = [];
%    end

   xi(i,:)=[3*randn(1) 2*randn(1)];     %generate 2 Principal components
                                        %1st PC score: N(0,9), 2nd PC score: N(0,4)
   %generate the repeated measurements with measurement errors
   y{i}=mu_true(t{i},lint)+xi(i,:)*xeig(t{i},lint,2)+randn(1,length(t{i}));
   newy{i} = mu_true(newt{i},lint) + xi(i,:)*xeig(newt{i},lint,2)+randn(1,length(newt{i}));
          
   %true mean function: t+sin(t)
   %1st eigenfunction: -sqrt(0.2)*cos(2*t*pi/10)
   %2nd eigenfuncton: sqrt0.2)*sin(2*t*pi/10)
   %measurement error is distributed as N(0,1)
end


%========================================================================================================
%New style of setting the input arguments for function FPCA() that calls PCA().
p = setOptions('yname','x', 'regular', regular, 'selection_k', 'FVE','FVE_threshold', 0.9,'screePlot',1, 'designPlot',1,'corrPlot',1,'numBins',0, 'verbose','on');  
%Here, the options "p" is set as the follows: yname = 'x', no. of PC is selected by FVE and FVE_threshold = 0.9,
%create scree plot and design plot (by default, they are not displayed), 
%do not perform prebinning for the input data (by default, the program try to bin the data)
%and the rest uses default values                           
%Note the default value of selection_k is BIC1, which is recommended.
%By default, we display the diagnostic messages. You can suppress the messages by setting 'verbose' to be 'off'.

%Use FPCA() to recover the functional object for y and the results is a cell array
%that can be assessible by getVal(), e.g., to get eigenfunctions "phi", use
%phi = getVal(yy, 'phi');
%This method is recommended!

fprintf(1,'Recover the individual trajectories using FPCA():\n');
time=cputime;
[yy] = FPCA(y,t,p);
time = cputime-time;
display(['Total time in FPCA fitting is : ' num2str(time) ' sec.']);

%extract some of the results for the plots below:
out1 = getVal(yy,'out1');      %vector of time points for mu, phi and ypred
mu = getVal(yy,'mu');          %estimated mean function
out21 = getVal(yy,'out21');    %vector of time points for xcov
xcovfit = getVal(yy,'xcovfit');%fitted covariance surface
xcov = getVal(yy,'xcov');      %estimated smooth covariance evaluated at out21
phi = getVal(yy,'phi');        %estimated eigenfunctions
no_opt = getVal(yy,'no_opt');  %number of optimal FPC scores
xcorr = getVal(yy,'xcorr');    %fitted correlation surface
%========================================================================================================

%create Kth mode variation plots for k = 1:no_opt
KModeVariationPlot(yy)  %kth mode variation plots for all k from 1 to no_opt
%Uncomment this if you just want k=1, use
%KModeVariationPlot(yy,1);

%Uncomment this to create scree plot after FPCA()
%createScreePlot(yy);

% compute the true covariance surface
truecov=[];
for i=1:length(out1)
    for j=1:length(out1)
        truecov(i,j)=xeig(out1(i),lint,2)'*diag([9,4])*xeig(out1(j),lint,2);
    end
end
ycov=truecov;
ycov = truecov + eye(length(out1));
truecorr=[];
for i=1:length(out1)
    for j=1:length(out1)
        truecorr(i,j)=truecov(i,j)/sqrt(truecov(i,i)*truecov(j,j));
    end
end
%for i=1:length(out1)
%    ycov(i,i)=truecov(i,i)+1;
%end
%===========Plot true covariance and diagonal variance===============================
figure;
set(gcf,'Position',[198 39 874 701]);
subplot(3,2,1)
plot(cell2mat(t),cell2mat(y),'k.');
hold on;
plot(out1,mu,'r');
plot(out1,mu_true(out1,lint))
xlabel('t');
ylabel('\mu(t)');
legend('data','fitted','true','Location','Best');
title('Mean function');
                                                                                                                                                             
subplot(3,2,2)
[a b]=meshgrid(out1,out1);
mesh(a,b,truecov);
axis([0 10 0 10 -2 3.5]);
hold on
mesh(a,b,ycov);
axis([0 10 0 10 -2 3.5]);
hold off
title('True covariance and diagonal variance');
xlabel('t');
ylabel('t');
subplot(3,2,3)
mesh(out21,out21,xcovfit);
xlabel('t');
ylabel('t');
title('Fitted covariance surface');
                                                                                                                                                             
subplot(3,2,4)
mesh(out21,out21,xcov);
xlabel('t');
ylabel('t');
title('Smooth covariance surface')
                                                                                                                                                             
subplot(3,2,5)
true_phi = xeig(out1,lint,2);
[phi_sign] = checkPhiSign(out1, phi(:,1), true_phi(1,:)');
plot(out1,phi_sign*phi(:,1),'r');
hold on
plot(out1,true_phi(1,:), 'b');
xlabel('t')
ylabel('\phi_1(t)')
legend('fitted','true','Location','Best');
title('\phi_1(t) Vs t');

subplot(3,2,6)
plot(out1,true_phi(2,:), 'b');
if no_opt > 1
    hold on
    [phi_sign] = checkPhiSign(out1, phi(:,2), true_phi(2,:)');
    plot(out1,phi_sign*phi(:,2),'r');
end
xlabel('t')
ylabel('\phi_2(t)')
title('\phi_2(t) Vs t');
if no_opt > 1
 legend('true','fitted','Location','Best');
else
 legend('true','Location','Best');
end
    

ask();

%===========Plot observed, true and predicted curves===================================
ypred=FPCAeval(yy,[],out1);         %obtain predicted curves for all existing subjects
                                    %also allows each subject evaluated at
                                    %different time points
%or simply use ypred = getVal(yy,'ypred'); %all subjects evaluated at the
%same out1

% For prediction of existing or new subjects, use  
% where yy have been obtained from FPCA(). newy and newt are cell arrays
% for new subjects, which have the same structure like y and t.  
[yprednew, xi_new, xi_var] =  FPCApred(yy, newy, newt);
%if all new subjects are evaluated at the same time, newt can be a row vector for 
%time points for one subject


xtrue=cell(1,ncohort);
for i=1:ncohort
    xtrue{i}=mu_true(out1,lint)+xi(i,:)*xeig(out1,lint,2);
end
%plot 4 randomly selected curves
figure
set(gcf,'Position',[1 29 1247 705]);
k = 1;

for i=mysample(1:ncohort,9,0)      %randomly sample 9 of the curves (without replacement) for plotting
    subplot(3,3,k)
    plot(t{i},y{i},'*',out1,xtrue{i},'b-',out1,ypred{i},'r-');
    title(['Subject ' num2str(i)]);
    xlabel('t');
    ylabel('X(t)');
    k= k+1;
    legend('obs','true','pred','Location','Best')
end
                                                                                                                                                             
