%  Add paths to data and functions

addpath ('..\..')

%  Last modified 26 March 2011

%  -----------------------------------------------------------------------
%                    Berkeley Growth Data
%  -----------------------------------------------------------------------

%  load growth.mat if it is available, otherwise, do the following
%  analyses

load growth

%  -----------------------------------------------------------------------

ncasem = 39;
ncasef = 54;
nage   = 31;

fid = fopen('hgtm.dat','rt');
hgtmmat = reshape(fscanf(fid,'%f'),[nage,ncasem]);

fid = fopen('hgtf.dat','rt');
hgtfmat = reshape(fscanf(fid,'%f'),[nage,ncasef]);

age = [ 1:0.25:2, 3:8, 8.5:0.5:18 ]';

%  -----------------------------------------------------------------------
%  Go to line 365 to load the results of analyses through monotone
%  smoothing, or to line 879 to load the results through registration
%  -----------------------------------------------------------------------

%  the dimensions of the data

ncasem = growth.ncasem;
ncasef = growth.ncasef;
nage   = growth.nage;

%  the heights of the 39 boys

hgtmmat = growth.hgtmmat;

%  the heights of the 54 girls

hgtfmat = growth.hgtfmat;

%  the 31 ages of measurement

age = growth.age;

%  --------------  Smooth the data non-monotonically  --------------
%  This smooth uses the usual smoothing methods to smooth the data,
%  but is not guaranteed to produce a monotone fit.  This may not
%  matter much for the estimate of the height function, but it can
%  have much more serious consequences for the velocity and
%  accelerations.  See the monotone smoothing method below for a
%  better solution, but one with a heavier calculation overhead.

%  -----------  Create fd objects   ----------------------------

%  A B-spline basis with knots at age values and order 6 is used

rng      = [1,18];
knots    = age;
norder   = 6;
nbasis   = length(knots) + norder - 2;
hgtbasis = create_bspline_basis(rng, nbasis, norder, knots);

%  --- Smooth these objects, penalizing the 4th derivative  --
%  This gives a smoother estimate of the acceleration functions

%  set up roughness penalty smoothing function smooth_basis

Lfdobj   = int2Lfd(4);
lambda   = 1e-1;
hgtfdPar = fdPar(hgtbasis, Lfdobj, lambda);

%  smooth the data

hgtmfd = smooth_basis(age, hgtmmat, hgtfdPar);
hgtffd = smooth_basis(age, hgtfmat, hgtfdPar);

%  plot the data and the smooth

plotfit_fd(hgtmmat, age, hgtmfd)
plotfit_fd(hgtfmat, age, hgtffd)

%  ---------------------------------------------------------------
%         Find minimum GCV value of lambda
%  ---------------------------------------------------------------

lnlam   = -6:0.25:0;
gcvsave = zeros(length(lnlam),1);
dfsave  = gcvsave;
for i=1:length(lnlam)
  hgtfdPari = fdPar(hgtbasis, Lfdobj, 10^lnlam(i));
  [hgtfdi, dfi, gcvi] = smooth_basis(age, hgtfmat, hgtfdPari);
  gcvsave(i) = sum(gcvi);
  dfsave(i)  = dfi;
end

%  plot the results

subplot(1,1,1)
phdl = plot(lnlam, gcvsave, 'k-o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}(\lambda)')
ylabel('\fontsize{13} GCV(\lambda)')

%  ---------------  plot the first 10 female data  ---------------

%  set up a fine mesh of ages for plotting

agefine = linspace(1,18,101)';

%  Height

hgtfmatfine = eval_fd(agefine, hgtffd(1:10));

phdl = plot(agefine, hgtfmatfine, '-');
set(phdl, 'LineWidth', 2)
hold on
plot(age, hgtfmat(:,1:10), 'o')
hold off
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height (cm)')
axis([1,18,60,200])

%  Velocity

velfmatfine = eval_fd(agefine, hgtffd(1:10), 1);

phdl = plot(agefine, velfmatfine, '-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height Velocity (cm/yr)')
axis([1,18,0,20])

%  Acceleration

accfmatfine = eval_fd(agefine, hgtffd(1:10), 2);

phdl = plot(agefine, accfmatfine, '-', ...
            [1,18], [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height Acceleration (cm/yr/yr)')
axis([1,18,-4,2])

%  height, velocity and acceleration

i=5;
subplot(3,1,1)
phdl = plot(agefine, hgtfmatfine(:,i), '-');
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} Height')
axis([1,18,60,200])
subplot(3,1,2)
phdl = plot(agefine, velfmatfine(:,i), '-');
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} Velocity')
axis([1,18,0,20])
subplot(3,1,3)
phdl = plot(agefine, accfmatfine(:,i), '-', ...
            [1,18], [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} Acceleration')
axis([1,18,-4,2])

%  plot velocities with knots at each age

plot((1:10),1)
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Velocity (cm/yr)')

%  plot velocities with 12 basis functions

hgtbasis12 = create_bspline_basis(rng, 12, norder);
hgtffd12   = smooth_basis(age, hgtfmat, hgtbasis12);

plot((1:10),1)
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Velocity (cm/yr)')

%  plot accelerations with 12 basis functions

plot((1:10),2)
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Acceleration (cm/yr^2)')

%  plot acceleration curves for the first 10 girls
%  estimated both by 12 basis
%  functions and by spline smoothing.

D2hgtfmat12 = eval_fd(agefine, hgtffd12(1:10), 2);

subplot(1,2,1)
plot(agefine, D2hgtfmat12, 'k-')
xlabel('\fontsize{12} Age')
ylabel('\fontsize{12} Acceleration (cm/yr^2)')
axis([1,18,-40,10])
axis('square')

subplot(1,2,2)
D2hgtfmat35 = eval_fd(agefine, hgtffd(1:10), 2);
plot(agefine, D2hgtfmat35, 'k-')
xlabel('\fontsize{12} Age')
ylabel('\fontsize{12} Acceleration (cm/yr^2)')
axis([1,18,-12,2])
axis('square')

print -dps2 'c:/MyFiles/fdabook1/figs.dir/twoaccelplots.ps'

subplot(1,1,1)
phdl = plot(agefine, D2hgtfmat35, 'k-', [1,18], [0,0], 'k:');
set(phdl, 'LineWidth', 1)
lhdl = line(agefine, mean(D2hgtfmat35,2));
set(lhdl, 'LineStyle', '--')
xlabel('\fontsize{13} Age')
ylabel('\fontsize{13} Acceleration (cm/yr^2)')
axis([1,18,-4,2])

%  ----------------------------------------------------------
%  Estimate standard error of measurement for velocity and 
%    acceleration, re-smooth using the reciprocal of variance
%    of estimate as a weight, and display results
%  ----------------------------------------------------------

%  set up function smooth_pos

norderse = 3;
nbasisse = nage + norderse - 2;
stderrbasis = create_bspline_basis([1,18], nbasisse, norderse, age);
Wfd0   = fd(zeros(nbasisse,1),stderrbasis);  %  initial value for Wfd

%  Males

hgtmfit      = eval_fd(age, hgtmfd);
hgtmres      = hgtmmat - hgtmfit;   %  residuals
hgtmresmnsqr = mean(hgtmres.^2,2);  %  mean squared residuals

%  positively smooth the mean squared residuals

Lfdobj = 1;             %  smooth toward a constant
lambda = 1e-3;          %  smoothing parameter
hgtfdPar = fdPar(Wfd0, Lfdobj, lambda);

Wfdm = smooth_pos(age, hgtmresmnsqr, hgtfdPar);

%  compute the variance and standard error of measurements

hgtmvar = eval_pos(age, Wfdm);
hgtmstd = sqrt(hgtmvar);

subplot(1,1,1)
plot(age, sqrt(hgtmresmnsqr), 'o', age, hgtmstd, 'b-')

%  update weight vector for smoothing data

wtvec = 1./hgtmvar;
wtvec = wtvec./mean(wtvec);

%  set up new smooth of the data using this weight vector

Lfdobj   = int2Lfd(4);
lambda   = 1e-2;
hgtfdpar = fdPar(hgtbasis, Lfdobj, lambda);

%  smooth the data again

[hgtmfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis(age, hgtmmat, hgtfdpar, wtvec);

%  display the results

growthdisplay(age, hgtmmat, hgtmfd, hgtmstd, y2cMap, 'male')

%  Females

hgtffit      = eval_fd(age, hgtffd);
hgtfres      = hgtfmat - hgtffit;   %  residuals
hgtfresmnsqr = mean(hgtfres.^2,2);  %  mean squared residuals

%  positively smooth the mean squared residuals

Wfdf = smooth_pos(age, hgtfresmnsqr, hgtfdPar);

%  compute the variance and standard error of measurements

hgtfvar = eval_pos(age, Wfdf);
hgtfstd = sqrt(hgtfvar);

subplot(1,1,1)
plot(age, hgtfresmnsqr, 'o', age, hgtfvar, 'b-')
plot(age, sqrt(hgtfresmnsqr), 'o', age, hgtfstd, 'b-')

plot(age, hgtfresmnsqr, 'ko', age, hgtfvar, 'k-')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Variance of Measurement')

print -dps2 'c:/MyFiles/fdabook1/figs.dir/growthvariance.ps'

%  update weight vector for smoothing data

wtvec = 1./hgtfvar;
wtvec = wtvec./mean(wtvec);

%  set up new smooth of the data using this weight vector

Lfdobj   = int2Lfd(4);
lambda   = 1e-1;
hgtfdPar = fdPar(hgtbasis, Lfdobj, lambda);

%  smooth the data again

[hgtffd, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis(age, hgtfmat, hgtfdPar, wtvec);

accffd = deriv(hgtffd,2);
accmat = eval_fd(agefine, accffd);
accmn  = mean(accmat(:,1:10),2);

plot(agefine, accmat(:,1:10), '-', [1,18], [0,0], 'r:')
lhdl = line(agefine, accmn);
set(lhdl, 'LineWidth', 2, 'LineStyle', '--', 'color', 'b')
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height Acceleration(cm/year/year)')
axis([1,18,-4,2])

%  display the results

growthdisplay(age, hgtfmat, hgtffd, hgtfstd, y2cMap, 'female')

%  ----------------------------------------------------------
%          Compute monotone smooths of the data  
%  ----------------------------------------------------------

%  ------  First set up a basis for monotone smooth   --------
%  We use b-spline basis functions of order 6
%  Knots are positioned at the ages of observation.

norder = 6;
nbasis = nage + norder - 2;
wbasis = create_bspline_basis(rng, nbasis, norder, knots);
Lfdobj = int2Lfd(3);  %  penalize curvature of velocity
lambda = 10^(-0.5);   %  smoothing parameter

% -----------------  Male data  --------------------

WfdParm = fdPar(fd(zeros(nbasis,ncasem),wbasis), Lfdobj, lambda);

[Wfdm, betam, hgtmhatfd] = smooth_monotone(age, hgtmmat, WfdParm);

% -----------------  Female data  --------------------

WfdParf = fdPar(fd(zeros(nbasis,ncasef),wbasis), Lfdobj, lambda);

[Wfdf, betaf, hgtfhatfd] = smooth_monotone(age, hgtfmat, WfdParf);

%  -----------  Save the analyses up to this point  ----------------

save growth_monotone_smooth

load growth_monotone_smooth

%  plot data and smooth, residuals, velocity, and acceleration

%  Males:

index = 1:ncasem;
for i = index
  beta      = betam(:,i);
  hgtmhat   = beta(1) + beta(2).*eval_mon(age, Wfdm(i));
  Dhgtmhat  =           beta(2).*eval_mon(age, Wfdm(i), 1);
  D2hgtmhat =           beta(2).*eval_mon(age, Wfdm(i), 2);
  subplot(2,2,1)
  plot(age, hgtmmat(:,i), 'go', age, hgtmhat, '-')
  axis([1, 18, 60, 200]);
  xlabel('Years');  title(['Height for male ',num2str(i)])
  resi = hgtmmat(:,i) - hgtmhat;
  subplot(2,2,2)
  plot(age, resi, '-o',     [1,18], [0,0], 'r--')
  axis([1,18,-1,1]);
  xlabel('Years');  title('Residuals')
  subplot(2,2,3)
  plot(age, Dhgtmhat, '-',  [1,18], [0,0], 'r--')
  axis([1,18,0,15]);
  xlabel('Years');  title('Velocity')
  subplot(2,2,4)
  plot(age, D2hgtmhat, '-')
  axis([1,18,-6,6]);
  xlabel('Years') ;  title('Acceleration')
  pause;
end

% Females:

index = 1:ncasef;
for i = index
  beta = betaf(:,i);
  hgtfhat   = beta(1) + beta(2).*eval_mon(age, Wfdf(i));
  Dhgtfhat  =           beta(2).*eval_mon(age, Wfdf(i), 1);
  D2hgtfhat =           beta(2).*eval_mon(age, Wfdf(i), 2);
  subplot(2,2,1)
  plot(age, hgtfmat(:,i), 'go', age, hgtfhat, '-')
  axis([1, 18, 60, 200]);
  xlabel('Years');  title(['Height for female ',num2str(i)])
  resi = hgtfmat(:,i) - hgtfhat;
  subplot(2,2,2)
  plot(age, resi, '-o',     [1,18], [0,0], 'r--')
  axis([1,18,-1,1]);
  xlabel('Years');  title('Residuals')
  subplot(2,2,3)
  plot(age, Dhgtfhat, '-')
  axis([1,18,0,15]);
  xlabel('Years');  title('Velocity')
  subplot(2,2,4)
  plot(age, D2hgtfhat, '-',  [1,18], [0,0], 'r--')
  axis([1,18,-6,6]);
  xlabel('Years') ;  title('Acceleration')
  pause;
end

%  Compute velocity functions over a fine mesh

velfmatfine = zeros(101,ncasef);
for i = 1:ncasef
  beta = betaf(:,i);
  velfmatfine(:,i) = beta(2).*eval_mon(agefine, Wfdf(i), 1);
end

velfmeanfine = mean(velfmatfine,2);

%  Compute acceleration functions over a fine mesh

accfmatfine = zeros(101,ncasef);
for i = 1:ncasef
  beta = betaf(:,i);
  accfmatfine(:,i) = beta(2).*eval_mon(agefine, Wfdf(i), 2);
end

accfmeanfine = mean(accfmatfine,2);

%  Plot the results for the first 10 girls

subplot(1,1,1)
index = 1:10;
accfmeanfine = mean(accfmatfine(:,index),2);
phdl = plot(agefine, accfmatfine(:,index), 'b-', ...
            [1,18], [0,0], 'b:');
set(phdl, 'LineWidth', 1)
lhdl = line(agefine, accfmeanfine);
set(lhdl, 'LineWidth', 2, 'LineStyle', '--')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Acceleration (cm/yr^2)')
axis([1,18,-4,2])

%  Phase-plane plot for the first 10 girls, 
%  Figure growloops in R-book

subplot(1,1,1)
index = 1:10;
phdl = plot(velfmatfine(:,index), -accfmatfine(:,index), 'b-', ...
            [1,18], [0,0], 'k:');
set(phdl, 'LineWidth', 1)
hold on
phdl = plot(velfmeanfine, -accfmeanfine, 'k--', ...
            [0,15], [0,0], 'k:');
set(phdl, 'LineWidth', 2)
phdl=plot(velfmatfine(62,index), -accfmatfine(62,index), 'ro', ...
          velfmeanfine(62), -accfmeanfine(62), 'ko');
set(phdl, 'LineWidth', 2)
text(1,7,'\fontsize{13} -- Mean curve')
text(1,6,'\fontsize{13} o','color','r')
text(1.5,6,'\fontsize{13} Age 11.3 (mean PGS ctr)')
hold off
xlabel('\fontsize{13} Velocity (cm/yr)')
ylabel('\fontsize{13} Deceleration (cm/yr^2)')
axis([0,15,-2,10])

%  ---------------------------------------------------------------------
%            Register the acceleration curves for the girls using
%                       landmark registration
%  ---------------------------------------------------------------------

load growth_registration

%  Define the mid-spurt as the single landmark by manually
%    clicking on the zero crossing

index = 1:ncasef;
midspurtsave = zeros(length(index),2);
D2hgtfmat    = zeros(ncasef, length(agefine));
subplot(1,1,1)
for i = index
  Wfdi  = fd(cvecf(:,i),wbasis);
  beta = betaf(:,i);
  D2hgtfmat(i,:) = beta(2).*eval_mon(agefine, Wfdi, 2);
  plot(agefine, D2hgtfmat(i,:), '-', [1,18], [0,0], 'r:')
  axis([1,18,-6,6]);
  xlabel('Years') ;  
  title(['Acceleration for record ',num2str(i)])
  midspurtsave(i,:) = ginput(1);
  pause;
end

%  This is an automatic PGS spurt identification procedure.
%  A mouse click advances the plot to the next case.
%  Compute PGS mid point for landmark registration.
%  Downward crossings are computed within the limits defined
%  by INDEX.  Each of the crossings within this interval 
%  are plotted.  The estimated PGS center is plotted as a vertical line.

%  The choice of range of argument values (6--18) to consider
%  for a potential mid PGS location is determined by previous
%  analyses, where they have a mean of about 12 and a s.d. of 1.

%  We compute landmarks for all 54 girls

accelfdUN = deriv(hgtffd, 2);

index  = 1:102;  %  wide limits
nindex = length(index);
ageval = linspace(8.5,15,nindex)';
PGSctr = zeros(ncasef,1);
subplot(1,1,1)
for icase = 1:ncasef
    accveci = eval_fd(ageval, accelfdUN(icase));
    aup     = accveci(2:nindex);
    adn     = accveci(1:(nindex-1));
    indx    = find(adn.*aup < 0 & adn > 0);
    plot(ageval(2:nindex),aup,'o',[8,18],[0,0],'b--')
    axis([7.9,18,-6,4])
    hold on
    for j = 1:length(indx)
        indxj = indx(j);
        aupj  = aup(indxj);
        adnj  = adn(indxj);
        agej  = ageval(indxj) + 0.1*(adnj./(adnj-aupj));
        if j == length(indx)
            PGSctr(icase) = agej;
            plot([agej,agej],[-4,4],'b-')
        else
            plot([agej,agej],[-4,4],'b:')
        end
    end
    hold off
    title(['Case ',num2str(icase)])
    pause
end

save PGSctr

PGSctrmean = mean(PGSctr);

subplot(1,1,1)
hist(PGSctr)
xlabel('\fontsize{16} Age of center of pubertal growth spurt')
ylabel('\fontsize{16} Frequency')
hold on
phdl=plot([PGSctrmean,PGSctrmean],[0,12],'r--');
set(phdl, 'LineWidth', 2)
hold off

%  compute the second derivative curve values over a fine mesh

%  Set up a simply monomial basis for landmark registration
%  This will compute warping functions that interpolate the
%  single landmark time

wbasisLM = create_bspline_basis([1,18], 4, 3, [1,PGSctrmean,18]);
WfdLM    = fd(zeros(4,1),wbasisLM);
WfdParLM = fdPar(WfdLM,1,1e-12);

%  carry out the landmark registration

[accelfdLM, warpfdLM, WfdLM] = ...
       landmarkreg(accelfdUN, PGSctr, PGSctrmean, WfdParLM, 1);
   
%  plot registered accelerations along with warping functions

accelmatUN = eval_fd(agefine, accelfdUN);
accelmatLM = eval_fd(agefine, accelfdLM);
warpmatLM  = eval_fd(agefine, warpfdLM);

for i=1:ncasef
    subplot(1,2,1)
    plot(agefine, accelmatLM(:,i), 'b-', ...
         agefine, accelmatUN(:,i), 'b--',...
         [1,18], [0,0], 'b:', ...
         [11.52,11.52], [-6,4], 'b:')
    axis([1,18,-6,4]);
    xlabel('Years')
    ylabel('Registered acceleration')
    axis('square')
    subplot(1,2,2)
    plot(agefine, warpmatLM(:,i), 'b-', ...
         [1,18], [1,18], 'b--', ...
         11.52, PGSctr(i), 'o')
    axis([1,18,1,18]);
    xlabel('Years')
    title(['Case ',num2str(i)])
    axis('square')
    pause
end

%  plot accelerations and warping functions for cases 3 and 7

m = 0;
for i=[3,7]
    m = m + 1;
    subplot(2,2,m)
    phdl=plot(agefine, accelmatUN(:,i), 'b-',...
         [1,18], [0,0], 'b:', ...
         [11.52,11.52], [-6,4], 'b--');
    set(phdl, 'LineWidth', 2)
    axis([1,18,-6,4]);
    m = m + 1;
    subplot(2,2,m)
    phdl=plot(agefine, warpmatLM(:,i), 'b-', ...
         [1,18], [1,18], 'b--', ...
         11.52, PGSctr(i), 'o', ...
         [11.5,11.5], [1,18], 'b--');
    set(phdl, 'LineWidth', 2)
    axis([1,18,1,18]);
end

%  plot registered accelerations for first 10 girls

figure(1)

%  plot the unregistered and registered results for the first 10 girls

subplot(2,1,1)
lhdl = plot(agefine, accelmatUN(:,1:10), '-', [1,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 1)
lhdl = line(agefine, mean(accelmatUN(:,1:10),2));
set(lhdl, 'LineWidth', 2, 'LineStyle', '--')
ylabel('Accel. (cm/yr^2)', 'FontSize', 13);
title('Unregistered', 'FontSize', 13);
axis([1,18,-4,2])

subplot(2,1,2)
lhdl = plot(agefine, accelmatLM(:,1:10), '-', [1,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 1)
lhdl = line(agefine, mean(accelmatLM(:,1:10),2));
set(lhdl, 'LineWidth', 2, 'LineStyle', '--')
xlabel('Age (years)', 'FontSize', 13);
ylabel('Accel. (cm/yr^2)', 'FontSize', 13);
title('Registered', 'FontSize', 13);
axis([1,18,-4,2])

%  plot the warping functions for the first 10 girls

subplot(1,1,1)
lhdl = plot(agefine, warpmatLM(:,1:10), '-', [1,18], [1,18], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age (clock time)')
ylabel('\fontsize{13} Age (growth time)')
axis([1,18,1,18])

%  plot the deformation functions for the first 10 girls

subplot(1,1,1)
lhdl = plot(agefine, warpmatLM(:,1:10)-agefine*ones(1,10), '-', ...
            [1,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age (clock time)')
ylabel('\fontsize{13} Age (growth retardation)')
axis([1,18,-3,3])

%  plot the warping and deformation functions for the first 10 girls

subplot(2,1,1)
lhdl = plot(agefine, warpmatLM(:,1:10), '-', [1,18], [1,18], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age (clock time)')
ylabel('\fontsize{13} Age (growth time)')
axis([1,18,1,18])

subplot(2,1,2)
lhdl = plot(agefine, warpmatLM(:,1:10)-agefine*ones(1,10), '-', ...
            [1,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age (clock time)')
ylabel('\fontsize{13} Retardation')
axis([1,18,-3,3])


%  ---------------------------------------------------------------------
%            Register the velocity curves for the girls using
%                     continuous registration
%  ---------------------------------------------------------------------

index = 1:ncasef;

%  set up a basis for the functions W(t) that define the warping
%  functions

nbasisw = 15;
norder  =  5;
basisw  = create_bspline_basis([1,18], nbasisw, norder);

%  define the target function for registration, as well as the
%  curves to be registered. 

Dhgtfmeanfd = deriv(mean(hgtffd), 1);
Dhgtffd     = deriv(hgtffd(index), 1);

%  define the functional parameter object for the W functions

coef0 = zeros(nbasisw,length(index));
Wfd0  = fd(coef0, basisw);

Lfdobj = int2Lfd(2);
lambda = 1;

WfdPar = fdPar(Wfd0, Lfdobj, lambda);

%  register the curves

tic;
[Dhgtfregfd, warpfd] = register_fd(Dhgtfmeanfd, Dhgtffd, WfdPar);
toc

tic;
[Dhgtfregfd, warpfd] = register_fd_QN(Dhgtfmeanfd, Dhgtffd, WfdPar);
toc

%  set up values for plotting

Dhgtfmat     = eval_fd(agefine, Dhgtffd);
Dhgtfregmat  = eval_fd(agefine, Dhgtfregfd);
Dhgtfmeanvec = eval_fd(agefine, Dhgtfmeanfd);
Dhgtfvec     = eval_fd(agefine, Dhgtffd);

%  set up values of the warping functions

warpmat = eval_fd(agefine, warpfd);

%  plot each curve

for i = 1:length(index)
   subplot(1,2,1)
   phdl = plot(agefine, Dhgtfvec(:,i),    '-',  ...
               agefine, Dhgtfmeanvec,     '--', ...
               agefine, Dhgtfregmat(:,i), '-',  ...
               [11.5, 11.5], [0,20], ':');
   set(phdl, 'LineWidth', 2)
   xlabel('\fontsize{13} Age (years)')
   ylabel('\fontsize{13} Growth velocity (cm/yr)')
   axis([1,18,0,20])
   axis('square')
   title(['\fontsize{13} Case ',num2str(i)])
   legend('\fontsize{12} Unregistered', ...
          '\fontsize{12} Target', ...
          '\fontsize{12} Registered')
   subplot(1,2,2)
   phdl = plot(agefine, warpmat(:,i), '-', agefine, agefine, '--');
   set(phdl, 'LineWidth', 2)
   xlabel('\fontsize{13} Growth age (years)')
   ylabel('\fontsize{13} Clock age (years)')
   axis([1,18,1,18])
   axis('square')
   pause
end

%  plot the unregistered and registered velocities for the first 10 girls

subplot(2,1,1)
lhdl = plot(agefine, Dhgtfmat(:,1:10), '-');
set(lhdl, 'LineWidth', 1)
lhdl = line(agefine, mean(Dhgtfmat(:,1:10),2));
set(lhdl, 'LineWidth', 2, 'LineStyle', '--')
ylabel('Velocity (cm/yr)', 'FontSize', 13);
title('Unregistered', 'FontSize', 13);
axis([1,18,0,20])

subplot(2,1,2)
lhdl = plot(agefine, Dhgtfregmat(:,1:10), '-');
set(lhdl, 'LineWidth', 1)
lhdl = line(agefine, mean(Dhgtfregmat(:,1:10),2));
set(lhdl, 'LineWidth', 2, 'LineStyle', '--')
xlabel('Age (years)', 'FontSize', 13);
ylabel('Velocity (cm/yr)', 'FontSize', 13);
title('Registered', 'FontSize', 13);
axis([1,18,0,20])

%  plot the warping functions for the first 10 girls

subplot(1,1,1)
lhdl = plot(agefine, warpmat(:,1:10), '-', [1,18], [1,18], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age (clock time)')
ylabel('\fontsize{13} Age (growth time)')
axis([1,18,1,18])

%  plot the deformation functions for the first 10 girls

subplot(1,1,1)
lhdl = plot(agefine, warpmat(:,1:10)-agefine*ones(1,10), '-', ...
            [1,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age (clock time)')
ylabel('\fontsize{13} Age (growth retardation)')
axis([1,18,-3,3])


%  ---------------------------------------------------------------------
%   Register the landmark registered  acceleration curves for the girls 
%                     using continuous registration
%  ---------------------------------------------------------------------

%  set up a basis for the functions W(t) that define the warping
%  functions

nbasisw = 15;
norder  =  5;
basisw  = create_bspline_basis([1,18], nbasisw, norder);

%  define the target function for registration, as well as the
%  curves to be registered. 

accelmeanfdLM = mean(accelfdLM);

%  define the functional parameter object for the W functions

coef0 = zeros(nbasisw,ncasef);
Wfd0  = fd(coef0, basisw);

Lfdobj = int2Lfd(1);
lambda = 1;
WfdPar = fdPar(Wfd0, Lfdobj, lambda);

%  register the curves

[accelfdCR, warpfdCR, WfdCR] = ...
                   registerfd(accelmeanfdLM, accelfdLM, WfdPar);

plot(accelfdCR)
axis([1,18,-4,3])

%  set up values for plotting

accelmatCR = eval_fd(agefine, accelfdCR);
warpmatCR  = eval_fd(agefine, warpfdCR);

%  plot each curve

for i = 1:length(index)
   subplot(1,2,1)
   phdl = plot(agefine, accelmatUN(:,i),    '-',  ...
               agefine, D2hgtfmeanvec,     '--', ...
               agefine, D2hgtfregmat(:,i), '-',  ...
               [11.5, 11.5], [-6,4], ':');
   set(phdl, 'LineWidth', 2)
   xlabel('\fontsize{13} Age (years)')
   ylabel('\fontsize{13} Growth acceleration (cm/yr^2)')
   axis([1,18,-6,4])
   axis('square')
   title(['\fontsize{13} Case ',num2str(i)])
   legend('\fontsize{12} Unregistered', ...
          '\fontsize{12} Target', ...
          '\fontsize{12} Registered')
   subplot(1,2,2)
   phdl = plot(agefine, warpmat(:,i), '-', agefine, agefine, '--');
   set(phdl, 'LineWidth', 2)
   xlabel('\fontsize{13} Growth age (years)')
   ylabel('\fontsize{13} Clock age (years)')
   axis([1,18,1,18])
   axis('square')
   pause
end

subplot(1,1,1)
lhdl = plot(agefine, Dhgtfregmat(:,1:10), 'k-', [1,18], [0,0], 'k:');
set(lhdl, 'LineWidth', 1)
lhdl = line(agefine, mean(Dhgtfregmat(:,1:10),2));
set(lhdl, 'LineWidth', 2, 'LineStyle', '--', 'color', 'k')
xlabel('Age (years)', 'FontSize', 13);
ylabel('Accel. (cm/yr^2)', 'FontSize', 13);
axis([1,18,0, 20])

save growth_registration

load growth_registration

save growth

%  ---------------------------------------------------------------------
%        Monotone smooth of short term height measurements
%  ---------------------------------------------------------------------

%  ---------------- input the data  ----------------------------------

clear;

fid  = fopen('onechild.dat','rt');
temp = fscanf(fid,'%f');
n    = 83;
data = reshape(temp, [n, 2]);
day  = data(:,1);
hgt  = data(:,2);
rng  = [day(1), day(n)];
wgt  = ones(n,1);
zmat = wgt;

%  set up the basis

nbasis   = 43;
norder   = 4;
hgtbasis = create_bspline_basis(rng, nbasis, norder);

%  set parameters for the monotone smooth

Lfdobj   = int2Lfd(2);
lambda   = 1;
hgtfdPar = fdPar(hgtbasis, Lfdobj, lambda);

%  carry out the monotone smooth

[Wfd, beta, Fstr, iternum, iterhist] = ...
                        smooth_monotone(day, hgt, hgtfdPar);

%  plot the function W = log Dh

subplot(1,1,1)
plot(Wfd);

dayfine  = linspace(day(1),day(n),151)';
yhat     = beta(1) + beta(2).*eval_mon(day, Wfd);
yhatfine = beta(1) + beta(2).*eval_mon(dayfine, Wfd);

%  plot the data plus smooth

phdl = plot(day, hgt, 'o', dayfine, yhatfine, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{19} Day')
ylabel('\fontsize{19} Height (cm)')
axis([0,312,123,131])

%  plot growth velocity

Dhgt = beta(2).*eval_mon(dayfine, Wfd, 1);
phdl = plot(dayfine, Dhgt);
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{19} Days')
ylabel('\fontsize{19} Velocity (cm/day)')
axis([0,312,0,.06])

%  plot growth acceleration

D2hgt = beta(2).*eval_mon(dayfine, Wfd, 2);
phdl = plot(dayfine, -D2hgt, [0,312], [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{19} Days')
ylabel('\fontsize{19} Deceleration (cm/day^2)')
axis([0,312,-.004,.003])

%  phase-plane plot

phdl = plot(Dhgt, -D2hgt, 'b-', [0,0.06], [0,0], 'r:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{19} Velocity (cm/day)')
ylabel('\fontsize{19} Deceleration (cm/day^2)')
axis([0,.06,-.004,.003])

