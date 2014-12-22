
%  Example file for how to use FPCder.m
%  Random design, 200 subjects with 1 to 8 measurements on each subject.
%  Time interval is [0,1]. The number of measurements for each subject is
%  uniformly distributed on {1,2,...,8}. The timepoints for each subject
%  are distributed as Beta(0.4,0.3).
%
%  In this example, the goal is to predict trajectories and their first two derivatives.

%generate data set
%clear all;

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
rand('twister',sum(10000*clock));
mtp = 8;        %at most 8 repeated measurements in the simulated data
ncohort=200;    %200 subjects in the simulated data
lint=1;
y=cell(1,ncohort);
t=cell(1,ncohort);
newy = y;
newt = t;
xi=zeros(ncohort,2);
ngrid = 100;

%Case iii) regular data with missing values (regular = 1)   
%ni_missing = zeros(1,ncohort);
%for i = 1:ncohort
%  ni_missing(i) = poissrnd(mysample(1:mtp,1,0),1,1); 
%  if ni_missing(i) >= mtp
%     ni_missing(i) = mtp-1;
% end
%end

for i=1:ncohort

   %Case i) Sparse and irregular case (regular = 0) 
   ntp=ceil(mtp*rand(1));
   t{i}=lint*rand(1,ntp);                 
   %newt{i} = lint*rand(1,ntp);
   newt{i} = sort(betarnd(0.4,0.3,1,ntp));

   %Case ii) complete balance case (regular = 2) 
   %t{i} = linspace(0,lint,mtp);                      
   %newt{i} = lint*rand(1,ntp);

   %Case iii) regular data with missing values (regular = 1)   
   %t{i} = linspace(0,lint,mtp);                    
   %newt{i} = linspace(0,lint,mtp);
   %if ni_missing(i) > 0
   %  id_missing = mysample(1:mtp, ni_missing(i),0);
   %  t{i}(id_missing) = [];
   %  newt{i}(mysample(1:mtp,ni_missing(i),0)) = [];
   %end

   xi(i,:)=[3*randn(1) 2*randn(1)];     %generate 2 Principal components
                                        %1st PC score: N(0,9), 2nd PC score: N(0,4)
   %generate the repeated measurements with measurement errors
   
   y{i}=genMeanFun_1(t{i},0)+xi(i,:)*genEigenFun_1(t{i},0)+randn(1,length(t{i}));
   newy{i} = genMeanFun_1(newt{i},0) + xi(i,:)*genEigenFun_1(newt{i},0)+randn(1,length(newt{i}));
          
   %measurement error is distributed as N(0,1)
end


 
%========================================================================================================

p = setDerOptions('yname','x','selection_k', 'FVE','FVE_threshold', 0.85,'screePlot',1, 'corrPlot',1,...
		  'designPlot',1,'numBins',0, 'nder',0:2,'newdata',linspace(0,1,ngrid), 'ngrid',ngrid,'verbose','on');  

% Here, the options "p" is set as the follows: yname = 'x', use FVE method for 
% selection of No. of FPCs (default is FVE, but if only the curve itself is of interest, set it as 'BIC1'), 
% FVE_threshold = 0.85, 
% create scree plot and design plot (by default, they are not displayed),
% do not perform prebinning for the input data (by default, the program try to bin the data)
% input data, perform curve estimation as well as 1st and 2nd derivative estimation 
% (by default nder = 0, only the curve itself will be estimated),
% define output time grids as 100 equidistant grids in [0,1],
% use 100 grids for the covariance surface and their partial derivative of 
% covariance surface, display diagnostic messages,
% and the rest uses default values                           
% It is important that 'nder' starts with 0, meaning the estimation of the trajectories themselves.
% By default, 'verbose' is set to be 'on'. You can suppress the diagnostic messages by setting it to be 'off'.

%Use FPCder() to recover the functional object for y and the results is a cell array
%that can be assessible by getVal(), e.g., to get eigenfunctions "phi", use
%phi = getVal(yy, 'phi');

fprintf(1,'Recover the individual trajectories using FPCder():\n');
time=cputime;
[yy] = FPCder(y,t,p);
time = cputime-time;
display(['Total time in FPCder fitting is : ' num2str(time) ' sec.']);


%create Kth mode variation plots for k = 1:no_opt
 
KModeVariationPlot(yy)  %kth mode variation plots for all k from 1 to no_opt

%Uncomment this if you just want k=1, use
%KModeVariationPlot(yy,1);
%Uncomment this to create scree plot after FPCder()
%createScreePlot(yy);

%extract some of the results for the plots below:

out1 = getVal(yy,'out1');      %vector of time points for mu, phi and ypred
mu = getVal(yy,'mu');          %estimated mean functions
mu0 = mu{1};                   %estimated mean function
mu1 = mu{2};                   %estimated 1st derivative of mean function
mu2 = mu{3};                   %estimated 2nd derivative of mean function

out21 = getVal(yy,'out21');    %vector of time points for xcov
xcovfit = getVal(yy,'xcovfit');%fitted covariance surface

xcov = getVal(yy,'xcov');         %estimated smooth covariances evaluated at out21
xcov0 = xcov{1};                  %estimated smooth covariances G(t,s)
xcov1 = xcov{2};                  %estimated 1st partial derivative of covariance G'(t,s)
xcov2 = xcov{3};                  %estimated 2nd partial derivative of covariance G''(t,s)


%Uncomment this to obtain the graph for fitted and true mean 
%functions and their first two derivatives
%compute the true mean function
%tmu0 = genMeanFun_1(out1,0);
%tmu1 = genMeanFun_1(out1,1);
%tmu2 = genMeanFun_1(out1,2);
%figure;
%subplot(3,1,1)
%plot(out1,mu0,'b-',out1,tmu0,'r-'); 
%xlabel('\bf{t}');
%ylabel('\bf{mu(t)}');
%legend('pred','true','Location','Best');
%subplot(3,1,2)
%plot(out1,mu1,'b-',out1,tmu1,'r-'); 
%xlabel('\bf{t}');
%ylabel('\bf{mu''(t)}');
%legend('pred','true','Location','Best');
%subplot(3,1,3)
%plot(out1,mu2,'b-',out1,tmu2,'r-'); 
%xlabel('\bf{t}');
%ylabel('\bf{mu''''(t)}');
%legend('pred','true','Location','Best');

% compute the true covariance surface
tcov0=genEigenFun_1(out21,0)'*diag([9,4])*genEigenFun_1(out21,0);  
ycov=tcov0;
ycov = tcov0 + eye(length(out21));  %add in variance due to measurement error to the diagonal
% compute the true partial derivative of covariance surface
%true 1st partial derivative of the covariance function 
tcov1 = genEigenFun_1(out21,0)'*diag([9,4])*genEigenFun_1(out21,1);
%true 2nd partial derivative of the covariance function  
tcov2 = genEigenFun_1(out21,0)'*diag([9,4])*genEigenFun_1(out21,2);  

%Uncomment this to obtain the graph for estimated smooth covariance
%function and its partial derivatives
%subplot(3,2,1)
%mesh(out21,out21,xcov0);
%xlabel('\bf{t}');
%ylabel('\bf{s}');
%zlabel('\bf{Estimated G(t,s)}');
%subplot(3,2,2)
%mesh(out21,out21,tcov0);
%xlabel('\bf{t}');
%ylabel('\bf{s}');
%zlabel('\bf{True G(t,s)}');
%subplot(3,2,3)
%mesh(out21,out21,xcov1);
%xlabel('\bf{t}');             %derivative direction
%ylabel('\bf{s}');             %smoothing direction
%zlabel('\bf{Estimated G''(t,s)}');
%subplot(3,2,4)
%mesh(out21,out21,tcov1);
%xlabel('\bf{t}');             %derivative direction
%ylabel('\bf{s}');             %smoothing direction
%zlabel('\bf{True G''(t,s)}');
%subplot(3,2,5)
%mesh(out21,out21,xcov2);
%xlabel('\bf{t}');             %derivative direction
%ylabel('\bf{s}');             %smoothing direction
%zlabel('\bf{Estimated G''''(t,s)}');
%subplot(3,2,6)
%mesh(out21,out21,tcov2);
%xlabel('\bf{t}');             %derivative direction
%ylabel('\bf{s}');             %smoothing direction
%zlabel('\bf{True G''''(t,s)}');



phi = getVal(yy,'phi');        %estimated eigenfunctions
phi0 = phi{1};                 %estimated eigenfunctions
phi1 = phi{2};                 %estimated 1st derivative of eigenfunctions
phi2 = phi{3};                 %estimated 2nd derivative of eigenfunctions
    
no_opt = getVal(yy,'no_opt');  %number of optimal FPC scores
%========================================================================================================


%===========Plot true covariance and diagonal variance===============================
figure;
set(gcf,'Position',[198 39 874 701]);
subplot(3,2,1)
plot(cell2mat(t),cell2mat(y),'k.');
hold on;
plot(out1,mu{1},'r');
plot(out1,genMeanFun_1(out1,0))
xlabel('t');
ylabel('\mu(t)');
legend('data','fitted','true','Location','Best');
title('Mean function');
      
subplot(3,2,2)
mesh(out21,out21,tcov0);
axis([0 10 0 10 -2 3.5]);
hold on
mesh(out21,out21,ycov);
axis([0 1 0 1 -20 20]);
hold off
title('True covariance and diagonal variance');
xlabel('t');
ylabel('t');
subplot(3,2,3)
mesh(out21,out21,xcovfit{1});
xlabel('t');
ylabel('t');
title('Fitted covariance surface');

subplot(3,2,4)
mesh(out21,out21,xcov{1});
xlabel('t');
ylabel('t');
title('Smooth covariance surface')
 
subplot(3,2,5)
true_phi = genEigenFun_1(out1,0);
[phi_sign] = checkPhiSign(out1, phi0(:,1), true_phi(1,:)');
plot(out1,phi_sign*phi0(:,1),'r');
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
    [phi_sign] = checkPhiSign(out1, phi0(:,2), true_phi(2,:)');
    plot(out1,phi_sign*phi0(:,2),'r');
end
xlabel('t')
ylabel('\phi_2(t)')
title('\phi_2(t) Vs t');
if no_opt > 1
 legend('true','fitted','Location','Best');
else
 legend('true','Location','Best');
end
hold off;                           
                                                                                                                                  
ask();

%===========Plot observed, true and predicted curves===================================
% For prediction of existing or new subjects, use
% where yy have been obtained earlier. newy and newt are cell arrays
% for new subjects, which have the same structure like y and t.
% ypred = FPCderPred(yy, y,out1,2); 
ypred = getVal(yy,'y_pred');
y0 = ypred{1};     %estimated curve
y1 = ypred{2};     %estimated first derivative
y2 = ypred{3};     %estimated second derivative

% Uncomment this if you are interested in prediction for new subjects
%[yprednew, xi_new, xi_var] =  FPCderPred(yy, newy, newt, []);


xtrue = repmat(genMeanFun_1(out1,0),ncohort,1)+xi*genEigenFun_1(out1,0);    %matrix of true x(t)
xtrue1 = repmat(genMeanFun_1(out1,1),ncohort,1)+xi*genEigenFun_1(out1,1); %matrix of true x'(t)
xtrue2 = repmat(genMeanFun_1(out1,2),ncohort,1)+xi*genEigenFun_1(out1,2); %matrix of true x''(t)


 ii=mysample(1:ncohort,4,0);      %randomly sample 4 of the curves (without replacement) for plotting
 for j = 1:ceil(length(ii)/2)
       figure
       set(gcf,'Position',[1 29 1247 705]);
       i = ii(j*2-1); 
  
       %estimated X(t)
       subplot(3,2,1);
       plot(t{i},y{i},'*',out1,xtrue(i,:),'b-',out1,y0{i},'r-');
       title(['Subject ' num2str(i)]);
       xlabel('t');
       ylabel('X(t)');
       legend('obs','true','pred','Location','Best');

       %estimated X'(t)
       subplot(3,2,3);
       plot(out1,xtrue1(i,:),'b-',out1,y1{i},'r-');
       title(['Subject ' num2str(i)]);
       xlabel('t');
       ylabel('X''(t)');
       legend('true','pred','Location','Best');

       %estimated X"(t)
       subplot(3,2,5);
       plot(out1,xtrue2(i,:),'b-',out1,y2{i},'r-');
       title(['Subject ' num2str(i)]);
       xlabel('t');
       ylabel('X''''(t)');
       legend('true','pred','Location','Best');

       i = ii(j*2);
       %estimated X(t)
       subplot(3,2,2);
       plot(t{i},y{i},'*',out1,xtrue(i,:),'b-',out1,y0{i},'r-');
       title(['Subject ' num2str(i)]);
       xlabel('t');
       ylabel('X(t)');
       legend('obs','true','pred','Location','Best');

       %estimated X'(t)
       subplot(3,2,4);
       plot(out1,xtrue1(i,:),'b-',out1,y1{i},'r-');
       title(['Subject ' num2str(i)]);
       xlabel('t');
       ylabel('X''(t)');
       legend('true','pred','Location','Best');

       %estimated X"(t)
       subplot(3,2,6);
       plot(out1,xtrue2(i,:),'b-',out1,y2{i},'r-');
       title(['Subject ' num2str(i)]);
       xlabel('t');
       ylabel('X''''(t)');
       legend('true','pred','Location','Best');
    
end

%plot dynamic transfer function \beta_1(t)
figure;
beta1 = getDynamic(yy);
plot([min(out1) max(out1)], [0 0], 'Color',[0.8 0.8 0.8]);
hold on;
plot(out1, beta1, 'b-');
xlabel('t');
ylabel('\beta_1(t)');
title('Dynamic transfer function \beta_1(t)');
hold off;
