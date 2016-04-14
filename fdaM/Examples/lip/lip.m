%%  Introduction to the lip position data with 20 replications

% These data record the position of the lower lip as a single speaker
% was saying "Bob" for each of 20 replications.  Also recorded at the same
% was the EMG (electromygraphical) activity in the principal muscle 
% controlling the verticalposition of the lower lip during speech. 
%
% The data were described and analyzed in Ramsay, Munhall, Gracco and Ostry 
% (1996), Functional data analyses of lip motion. Journal of the Acoustical 
% Society of America, 99, 3718-3727.  
%
% We were interested in exploring the dyanamics of lip movement, and 
% especially its acceleration in response to muscle activity.
% Some of our analyses looked only at the response, and others studied the
% relationship between the input EMG signal and output position or 
% acceleration.

%%
% Time was recorded in milliseconds, and the data in millimetres.
% But the choices of scales for both the observations and for time require
% some comment on computational issues and other comments on the best
% scales for interpreting the results.
% For example, can be tempting to change the time scale to seconds, so that
% the width of the observational interval is 0.35.   But the effect on the
% analyses due to rounding error is catastrophic.  You can try T - 0.35
% below to see for yourself what happens.
% It's best in general to choose a time scale such that the time difference
% between successive knots for B-spline basis functions does not go much 
% below 1.  For time in seconds, this difference becomes 0.35/50 = 0.007, 
% and has a catastrophic effect on estimation and analysis of acceleration.  
% Expressing time in 100th's of a second or in milliseconds is better.
%
% The observational scale of millimetres works well for plotting the 
% curves.  But it comes to plotting acceleration, it's better to
% plot on a scale that better related to our everyday experience.  The
% acceleration due to gravity is about 10 metres/second/second, and the
% various twitches and other sharp movements of the body, whether
% intentional or otherwise, seldom exceed about 30 metres/second/second.
% Our bodies are not engineered to survive much more acceleration than
% that.  Switching to this scale by changing the observational units
% makes for both better plotting and a more interpretable display.
%
% In these analyses, we use B-spline basis functions, which are the usual
% choice when the data are not periodic.  Since we want to study lip 
% acceleration, we want to control the smoothness of the second derivative,
% and this argues for the use of order 6 B-splines, rather than the 
% common choice of order 4, useful when only the smoothness of the curve
% itself is of interest.

%  Last modified 22 December 2014

%% Setting up the data for analysis

addpath('../..')  % path to location of functional data analysis software

fid    = fopen('lip.dat','rt');
lipmat = reshape(fscanf(fid,'%f'), [51, 20]);

%% Time in 100th's of a second
% Set up the b-spline basis object, using 56 order 6 splines so that there
% is a knot located at each time of observation.  The number of knots is
% the number of interior points (49) plus the norder (6), or 55.

[n,N] = size(lipmat);  % n = no. obs., N = no. records

T = 35;
liptime  = (0:n-1)'*T/(n-1); 

%  set up a fine time scale for plotting purposes

nfine = 500;
liptimefine = (0:nfine)*T/nfine;

%  set up the B-spline basis system

liprng = [0, T];  %  normalized time interval in seconds
lipnbasis = 55;
lipnorder =  6;
lipbasis  = create_bspline_basis(liprng, lipnbasis, lipnorder);

%  set up the functional parameter object

% penalize fourth derivative for smooth acceleration estimate
Lfdobj = int2Lfd(4);  
% Although the data are highly accurate, use fairly heavy smoothing  
% in order to estimate a smooth second derivative.
lambda = 1e2;      
% define the functional parameter object
lipfdPar = fdPar(lipbasis, Lfdobj, lambda);

%  carry out the smoothing to construct the functional data object lipfd
 
lipfd = smooth_basis(liptime, lipmat, lipfdPar);

%  add names to dimensions

lipfd_fdnames{1} = 'time (seconds*100)';
lipfd_fdnames{2} = 'Replications';
lipfd_fdnames{3} = 'mm';
lipfd = putnames(lipfd, lipfd_fdnames);

%  plot data and fit for each curve.  This cycles through each curve,
%  and is commented out to make it optional

% plotfit_fd(lipmat, liptime, lipfd)

%%  ------  summarize and plot the functions and their accelerations  -----

%  compute the mean functional data object

lipmeanfd = mean(lipfd);

%  evaluate lipfd and its second derivative

lipmat   = eval_fd(liptime, lipfd);
D2lipmat = eval_fd(liptime, lipfd, 2);

%  plot the evaluated values

figure(1)
subplot(2,1,1)
phdl = plot(liptime, lipmat, '-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Time (sec*100)')
ylabel('\fontsize{13} Lip position (mm)')

%  convert to metres per second per second by multiplying the data by 10

subplot(2,1,2)
phdl = plot(liptime, D2lipmat*10, '-', liprng, [0,0], 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Time (sec*100)')
ylabel('\fontsize{13} Lip acceleration (metres/sec^2)')

%%  Register the data 

%  Examining the accleration curves, we see that there is a certain amount
%  of variation across the curves in prominent features.  We call this
%  phase variation.

%  We now want to warp time so that curve features line up across curves.
%  We can easily identify two landmarks, the two maxima in each curve.
%  Our first registration procedure is to line up these two landmarks.

nmarks = 2;

% if the landmarks have already been identified and saved as below,
%  skip over the following to the "load landmarks" statement

%  set up the matrix of acceleration values for identifying landmarks

%  plot each acceleration curve, and click on the two
%    maxima, near t = 0.12 seconds, the other near t = .27 seconds

%  locate landmarks manually by clicking the minimum and elbow
%  for each curve

lipmarks = zeros(nobs,nmarks);
subplot(1,1,1)
for i = 1:nobs
    phdl = plot(liptime, D2lipmat(:,i)/1000, 'b-', liprng, [0,0], 'b:');
    axis([liprng,-15,10])
    set(phdl, 'LineWidth', 2)
    title(['\fontsize{16} Record ',num2str(i)])
    hold on
    [x,y] = ginput(1);
    for j = 1:nmarks
        [x,y] = ginput(1);  %  input clicks on points here
        lipmarks(i,j) = x(1);
    end
    hold off
    pause
end

% save lip marks for future work

save lipmarks lipmarks 

%  load the landmarks if they have already been saved

load lipmarks 

lipmeanmarks = mean(lipmarks);  %  mean lipmarks

%  plot each record along with the location of the two maxima

figure(2)
subplot(1,1,1)
for i=1:20
    phdl = plot(liptime, D2lipmat(:,i)*10, '-', liprng, [0,0], 'b--');
    set(phdl, 'LineWidth', 2)
    hold on
    for k=1:2
        phdl = plot([lipmarks(i,k),lipmarks(i,k)], [-15,10], 'b--');
        set(phdl, 'LineWidth', 2)
    end
    hold off
    xlabel('\fontsize{16} Time (sec*100)')
    ylabel('\fontsize{16} Lip acceleration (metres/sec^2)')
    title(['\fontsize{16} Record ',num2str(i)])
    pause
end

%%  Landmark registration of the curves  
%  In the first registration, the monotone smoothing option 
%  is used to compute the warping function.  This is slower,
%  but gives better results because the warping functions
%  are smoother.  See  below for a smooth not using the 
%  monotone smoothing option.

%  create a basis object for the warping function
%  it has order 4 (piecewise cubic) and two interior knots
%  positioned at the mean landmark values since
%  NBASIS = NORDER + # interior knots

nbasis = 13;
norder =  4;
breaks = [0,lipmeanmarks,T];
warpbasis = create_bspline_basis(liprng, nbasis, norder, breaks);
WfdPar    = fdPar(warpbasis,int2Lfd(2),1e-8);

%  plot the basis

subplot(1,1,1)
plot(warpbasis)  %  plot of B-spline basis functions
line([lipmeanmarks(1),lipmeanmarks(1)],[0,1])  % first knot
line([lipmeanmarks(2),lipmeanmarks(2)],[0,1])  % second knot

%  call landmark registration function to set up struct LMRKSTR

[lipregfd, lipwarpfd, Wfd] = ...
     landmarkreg(lipfd, lipmarks, lipmeanmarks, WfdPar,0,1e-8);

%  plot unregistered and registered position curves

lipregmat   = eval_fd(liptime, lipregfd);
D2lipregmat = eval_fd(liptime, lipregfd, 2);

figure(3)
subplot(2,1,1)
plot(liptime, lipmat, '-', ...
     [lipmeanmarks(1),lipmeanmarks(1)], [-15,15], 'b--', ...
     [lipmeanmarks(2),lipmeanmarks(2)], [-15,15], 'b--')
ylabel('\fontsize{16} Position (mm)')
title('\fontsize{16} Unregistered')

subplot(2,1,2)
plot(liptime, lipregmat, '-', ...
     [lipmeanmarks(1),lipmeanmarks(1)], [-15,15], 'b--', ...
     [lipmeanmarks(2),lipmeanmarks(2)], [-15,15], 'b--')
xlabel('\fontsize{16} Time (secs*100)')
ylabel('\fontsize{16} Position (mm)')
title('\fontsize{16} Landmark Registered')

%  plot unregistered and registered acceleration curves

figure(4)
subplot(2,1,1)
plot(liptime, D2lipmat*10, '-', ...
     [lipmeanmarks(1),lipmeanmarks(1)], [-15,15], 'b--', ...
     [lipmeanmarks(2),lipmeanmarks(2)], [-15,15], 'b--')
axis([0,T,-15,15])
ylabel('\fontsize{16} Acceleration (metres/sec^2)')
title('\fontsize{16} Unregistered')

subplot(2,1,2)
plot(liptime, D2lipregmat*10, '-', ...
     [lipmeanmarks(1),lipmeanmarks(1)], [-15,15], 'b--', ...
     [lipmeanmarks(2),lipmeanmarks(2)], [-15,15], 'b--')
axis([0,T,-15,15])
xlabel('\fontsize{16} Time (secs*100)')
ylabel('\fontsize{16} Acceleration (metres/sec^2)')
title('\fontsize{16} Landmark Registered')

%  plot warping functions h(t) and deformation functions: d(t) = h(t) - t

lipwarpmat = eval_fd(liptime, lipwarpfd);
figure(5)
subplot(2,1,1)
phdl = plot(liptime, lipwarpmat, '-');
set(phdl, 'LineWidth', 2)
axis([0,T,0,T])
xlabel('\fontsize{16} Clock time (secs*100)')
ylabel('\fontsize{16} Warped time h(t) (secs*100)')
subplot(2,1,2)
lipwarpmat = eval_fd(lipwarpfd,liptime);
lipdefmat  = lipwarpmat - liptime*ones(1,N);
phdl = plot(liptime, lipdefmat, '-', liprng, [0,0], ':');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Time (secs*100)')
ylabel('\fontsize{16} h(t) - t (secs*100)')

%%  Principal component analysis of the position curves
%  We use six principal components,  see log10 eigenvalue plot below
%  for justification

nharm = 6;
lippcastr = pca_fd(lipfd, nharm);

%  plot unrotated harmonics, cycle through by hitting any key

figure(6)
subplot(1,1,1)
plot_pca_fd(lippcastr)

%  rotate harmonics using the Varimax criterion

lippcastr = varmx_pca(lippcastr);

%  plot rotated harmonics

figure(6)
plot_pca_fd(lippcastr)
  
%  Plot log eigenvalues, note that a break seems to appear between
%  harmonics 6 and 7, motivating our choice of 6 harmonics

figure(7)
lipeigvals = lippcastr.values;
phdl = plot(1:19,log10(lipeigvals(1:19)),'-o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Eigenvalue Number')
ylabel('\fontsize{16} Log10 Eigenvalue')

%% Principal differential analysis 
%  Here we estimate a second order linear differential equation whose 
%  two-dimensional solution space best approximates the observed curves.

%  Estimate the weight functions w_1(t) and w_2(t) 

difeorder = 2;  %  order of equation

%  external forcing function, in this case not present

awtcell = {};
ufdcell = {};

%  set up a functional parameter object for the weight functions

nwbasis = 21;
wbasis = create_bspline_basis(liprng, nwbasis);
wcoef0 = zeros(nwbasis,1);
betafd = fd(wcoef0, wbasis); 
betafdPar = fdPar(betafd);

%  for the weight function parameters and the observed curves in 
%  cell arrays

bwtcell{1,1} = betafdPar;
bwtcell{1,2} = betafdPar;

xfdcell{1} = lipfd;

%  carry out principal differential analysis

[bfdcell, resfdcell, afdcell] = ...
    pda_fd(xfdcell, bwtcell, awtcell, ufdcell, difeorder);

%  plot the two weight functions, labelled 0 and 1, respectively

figure(8)
for j=1:2
    subplot(2,1,j)
    bfdj = getfd(bfdcell{1,j});
    bvecj = eval_fd(liptime, bfdj);
    phdl = plot(liptime, bvecj, '-', [0,T],[0,0], 'b--');
    set(phdl, 'LineWidth', 2)
    xlabel('\fontsize{16} Time (secs*100)')
    ylabel(['\fontsize{16} Weight function ',num2str(j-1)]);
end

%%  Evaluate the fit to the observed curves offered by the 
%   solutions of the second ordfer differential equation.

%  set up the estimated linear differential operator lipLfd

wcoef   = [getcoef(getfd(bfdcell{1})), getcoef(getfd(bfdcell{2}))];
wfd     = fd(wcoef,wbasis);
bwtcell = fd2cell(wfd);
lipLfd  = Lfd(difeorder, bwtcell);

%  Compute forcing functions, these are effectively residuals from the 
%  fits to the observed curves

force = eval_fd(liptime, lipfd, lipLfd);

%  Plot the forcing functions for each curve.  These are the differential
%  analogues of residual functions in regression analysis.

figure(9)
subplot(1,1,1)
plot(liptime, force, '-', [0,T], [0,0], 'b--');
xlabel('\fontsize{16} Time (secs*100)')
ylabel('\fontsize{16} Forcing functions Lx(t)');

%  Plot the mean forcing function along with second derivative in order
%  see how small the second order forcing functions are with respect 
%  to the simple second derivative.

forcemean = mean(force,2);
D2mean = eval_fd(mean(lipregfd),liptime,int2Lfd(2));

figure(10)
subplot(1,1,1)
phdl = plot(liptime, forcemean, '-', liptime, D2mean, 'b:');
set(phdl, 'LineWidth', 2)
axis([0,T,-0.6,0.6])
xlabel('\fontsize{16} Time (secs*100)')
ylabel('\fontsize{16} Forcing function Lx(t)');
legend('\fontsize{16} Mean Lx(t)', '\fontsize{16} Mean D^2x(2)')

%%  Solve estimated differential equation for two sets of initial
%   values, and use these two solutions to fit each curve by
%   computing the best approximation based on the two solutions.

ystart = eye(2);
[tp1, yp1] = ode45(@derivs_pda, liptime, ystart(:,1), [], bwtcell);
[tp2, yp2] = ode45(@derivs_pda, liptime, ystart(:,2), [], bwtcell);

%  Plot the two sets of solutions

figure(11)
umat = [yp1(:,1),yp2(:,1)];
subplot(2,1,1)
phdl = plot(liptime, umat, liptime, zeros(51,1), ':'); 
title('\fontsize{16} Function');
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{16} First solutions');
Dumat = [yp1(:,2),yp2(:,2)];
subplot(2,1,2)
phdl = plot(liptime, Dumat, liptime, zeros(51,1), ':');
title('\fontsize{16} Derivative');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Time (secs*100)')
ylabel('\fontsize{16} Second solutions');

%  Plot fit to each curve along with the associated forcing function and
%  the second derivative function as a benchmark

figure(12)
index = 1:N;
lipmat   = eval_fd(liptime, lipfd);
D2lipmat = eval_fd(liptime, lipfd, 2);
forcemat = eval_fd(liptime, lipfd, lipLfd);
for i = index
    
   subplot(2,1,1)
   %  solid is fit, dashed is actual
   xhat = umat * (umat\lipmat(:,i));
   phdl = plot(liptime, xhat, 'b-', liptime, lipmat(:,i), 'bo');
   set(phdl, 'LineWidth', 2)
   xlabel('\fontsize{16} Time (secs*100)')
   ylabel('\fontsize{16} Function x(t)')
   legend('\fontsize{16} Fitted x(t)', '\fontsize{16} Data', ...
          'Location', 'North')
      
   subplot(2,1,2)
   %  solid line is forcing function, dashed 2nd deriv.
   phdl = plot(liptime, forcemat(:,i), 'b-',  ...
               liptime, D2lipmat(:,i), 'b:', ...
               liprng, [0,0],   'r:');
   set(phdl, 'LineWidth', 2)
   axis([0,35,-1,1])
   title(['\fontsize{16} Record ',num2str(i)])
   ylabel('\fontsize{16} Forcing function u(t)')
   legend('\fontsize{16} u(t)', '\fontsize{16} D^2x(t)', ...
          'Location', 'South')
      
   pause;
end

%  save all of the results to this point

save lip


