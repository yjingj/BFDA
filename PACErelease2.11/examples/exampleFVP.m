
%  Example file for how to use FVP.m.
%  Random design, 200 subjects with 21 measurements on each subject.
%  Time interval is [0,10]. The timepoints for each subject are 0:.5:10. 
%

%generate data set
%clear all;

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
rand('twister',sum(10000*clock));
mtp = 21;        %21 repeated measurements in the simulated data
ncohort=200;     %200 subjects in the simulated data
lint=10;
y=cell(1,ncohort);
t=cell(1,ncohort);
r=y;
xi=[];  rxi =[];


for i=1:ncohort

   %complete balance data (regular = 2) 
   t{i} = linspace(0,lint,mtp);                      

   xi(i,:)=[3*randn(1) 2*randn(1)];     %generate 2 Principal components
                                        %1st PC score: N(0,1), 2nd PC score: N(0,.25)
   rxi(i,:)=[.5*randn(1) .3*randn(1)];     %generate 2 Principal components
                                        %1st PC score: N(0,1), 2nd PC score: N(0,.25)
   %generate the repeated measurements with measurement errors
   r{i}=(.2*mu_true(t{i},lint)+rxi(i,:)*xeig(t{i},lint,2)+.2*randn(1,length(t{i})));
   y{i}=mu_true(t{i},lint)+xi(i,:)*xeig(t{i},lint,2)+sign(randn(1,mtp)).*exp(r{i}./2);
          
   %true mean function: t+sin(t)
   %1st eigenfunction: -sqrt(0.2)*cos(2*t*pi/10)
   %2nd eigenfuncton: sqrt0.2)*sin(2*t*pi/10)
   %measurement error is distributed as N(0,9)
end

p = setOptions('yname','x','regular',2,'selection_k', 'FVE','FVE_threshold', 0.9,'screePlot',0, 'designPlot',0, 'verbose','on');  
[rr,sigma,yy] = FVP(y,t,p);

%extract some of the results for the plots below:
out1 = getVal(rr,'out1');      %vector of time points for mu, phi and ypred
rpred=FPCAeval(rr,[],out1);

id = mysample(1:ncohort,9,0);
%===========Plot true and predicted FVP curves===================================
figure;
% set(gcf,'Position',[1 29 1247 705]);
k=1;
for i=id      %randomly sample 9 of the curves (without replacement) for plotting
    subplot(3,3,k)
    plot(t{i},r{i},'.',out1,rpred{i},'--');
    title(['Subject ' num2str(i)]);
    xlabel('t');
    ylabel('R(t)');
    k= k+1;
end
legend('True','FVP','Location','Best')

