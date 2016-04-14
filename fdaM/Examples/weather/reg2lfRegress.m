

addpath ('../..')

%  Last modified 22 April 2009

% load TempClimateZone

%  ---------------------------------------------------------------
%             Predicting temperature from climate zone 
%  ---------------------------------------------------------------

fid    = fopen('dailtemp.dat','rt');
tempav = fscanf(fid,'%f');
tempav = reshape(tempav, [365,35]);

smallnbasis = 65;
daytime   = (1:365)'-0.5;
dayrange  = [0,365];
dayperiod = 365;
smallbasis  = create_fourier_basis(dayrange, smallnbasis);
[tempfd, df, gcv, beta, SSE, penmat, y2cMap, argvals, y] = smooth_basis(daytime, tempav, smallbasis);

place = [ ...
    'Arvida     '; ...
    'Bagottville'; ...
    'Calgary    '; ...
    'Charlottown'; ...
    'Churchill  '; ...
    'Dawson     '; ...
    'Edmonton   '; ...
    'Fredericton'; ...
    'Halifax    '; ...
    'Inuvik     '; ...
    'Iqaluit    '; ...
    'Kamloops   '; ...
    'London     '; ...
    'Montreal   '; ...
    'Ottawa     '; ...
    'Pr. Albert '; ...
    'Pr. George '; ...
    'Pr. Rupert '; ...
    'Quebec     '; ...
    'Regina     '; ...
    'Resolute   '; ...
    'Scheffervll'; ...
    'Sherbrooke '; ...
    'St. Johns  '; ...
    'Sydney     '; ...
    'The Pas    '; ...
    'Thunder Bay'; ...
    'Toronto    '; ...
    'Uranium Cty'; ...
    'Vancouver  '; ...
    'Victoria   '; ...
    'Whitehorse '; ...
    'Winnipeg   '; ...
    'Yarmouth   '; ...
    'Yellowknife'];

%  names for climate zones

zonename = ['Atlantic'; 'Pacific '; 'Contintl'; 'Arctic  '];

%  indices for weather stations in each of four climate zones

for i=1:35
    display([num2str(i), ' ', place(i,:)])
end

atlindex = [1,2,4,8,9,13,14,15,19,23,24,25,28,34];
pacindex = [12,17,18,30,31];
conindex = [3,5,6,7,16,20,26,27,29,32,33,35];
artindex = [10,11,21,22];

%  Set up a design matrix having a column for the grand mean, and
%    a column for each climate zone effect. Add a dummy contraint
%    observation

zmat = zeros(35,4);
zmat(atlindex,1) = 1;
zmat(pacindex,2) = 1;
zmat(conindex,3) = 1;
zmat(artindex,4) = 1;

p = 4;
xfdcell = cell(1,p);
for j=1:p
    xfdcell{j} = zmat(:,j);
end

%  set up the basis for the regression functions
nbetabasis = 13;
betabasis  = create_fourier_basis(dayrange, nbetabasis);

%  ----  set up the harmonic acceleration operator  -------
Lbasis  = create_constant_basis(dayrange);  %  create a constant basis
Lcoef   = [0,(2*pi/365)^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object


%  set up the functional parameter object for the regression fns.

betafd    = fd(zeros(nbetabasis,p), betabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betafd, harmaccelLfd, lambda, estimate);
betacell = cell(1,p);
for j=1:p
    betacell{j} = betafdPar;
end

%  initialize iterations

basisw   = create_bspline_basis([0,365], 5, 4);
% basisw   = create_bspline_basis([0,365],1,1);
niter    = 2;
itercell = cell(niter,4);
nharm    = 1;
WfdPar   = fdPar(basisw);
periodic = 1;

N = 35;
xfd      = tempfd; % fd of raw data
xmeanfd  = mean(xfd);
xmeanmat = eval_fd(daytime, xmeanfd)*ones(1,N);

nfine = length(daytime);
tmat = daytime*ones(1,N);

%  carry out the iterative registration

for iter=1:niter
    
    display(['************  Iteration ',num2str(iter), ...
             ' ************'])
    
    %  regression analysis
    
    fRegressCell = fRegress(tempfd, xfdcell, betacell, ones(1, N), y2cMap);
    
    %  construct model matrix
    
    betaestcell = fRegressCell.betahat;
    betastderr_cell = fRegressCell.betastderr;
    xhatfd      = fRegressCell.yhat; % fited fd data
    
    %  register to the model
    
    [yfd, warpfd, Wfd, shift] = ...
          register_fd(xhatfd, xfd, WfdPar, periodic);
    
    %  compute deformations
    
    deformmat = eval_fd(daytime, warpfd) - tmat;
    
    %  save results
    
    itercell{iter, 1} = yfd;
    itercell{iter, 2} = deformmat;
    itercell{iter, 3} = betaestcell;
    itercell{iter, 4} = shift;
    
    %  plot results
    
    figure(1)

    subplot(2,1,1)
    plotfit_fd(tempav, daytime, xfd, [], [], [], 1)
    xlabel('')
    ylabel('\fontsize{13} xfd')
    subplot(2,1,2)
    plotfit_fd(tempav, daytime, yfd, [], [], [], 1)
    xlabel('')
    ylabel('\fontsize{13} yfd')
    
    figure(2)
    plot(daytime, deformmat, '-', [0,365], [0,0], 'r:')
    
    figure(3) 
    plot_beta(betaestcell{1})
    hold on
    for j=2:4
        plot_beta(betaestcell{j})
    end
    axis([0,365,-30,20])
    
    pause;
    
    %  set up for the next iteration
    xfd      = yfd;
    xmeanfd  = mean(xfd);
    xmeanmat = eval_fd(daytime, xmeanfd)*ones(1,N);
    
end

%  plot final results

%  xfd and yfd's

figure(1)
subplot(1,1,1)
xmat = eval_fd(daytime, tempfd);
yarray = zeros(nfine,N,niter);
for iter=1:niter
    yarray(:,:,iter) = eval_fd(daytime, itercell{iter,1});
end
for i=1:N
    effecti = zmat(i,:);
    zonei   = find(effecti == 1);
    betafdPari = betaestcell{zonei};
    betafdi    = getfd(betafdPari);
    betaveci   = eval_fd(daytime, betafdi);
    plot(daytime, xmat(:,i), 'g-', daytime, betaveci, 'r--', ...
         [0,365], [0,0], 'r:')
    hold on
    for iter=1:niter
        plot(daytime, yarray(:,i,iter), '-')
    end
    hold off
    axis([0,365,-35,25])
    xlabel('\fontsize{13} Age (years)')
    ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
    title(['\fontsize{16} Case ',num2str(i),' ',place(i,:)])
    pause
end

%  regression functions

figure(2)
for j=1:p
    betacell1 = itercell{1,3};
    betafdj1 = getfd(betacell1{j});
    plot(betafdj1)
    hold on
    for iter=2:niter
        betacelli = itercell{iter,3};
        betafdji  = getfd(betacelli{j});
        plot(betafdji)
    end
    hold off
    axis([0,365,-35,25])
    xlabel('\fontsize{13} Age (years)')
    ylabel('\fontsize{13} Regression function')
    title(['\fontsize{16} ',zonename(j,:)])
    pause
end

%  deformation functions

deformtotal = zeros(nfine, N);
shifttotal  = zeros(N,1);
for iter=1:niter
    deformtotal = deformtotal + itercell{iter,2};
    shifttotal  = shifttotal  + itercell{iter,4};
end

display([min(shifttotal),max(shifttotal)])

figure(3)
for i=1:N
    plot(daytime, deformtotal(:,i), 'g-', [0,365], [0,0], 'r:')
    hold on
    for iter=1:niter
        plot(daytime, itercell{iter,2}(:,i), '-')
    end
    hold off
    xlabel('\fontsize{13} Age (years)')
    ylabel('\fontsize{13} Deformation h(t) - t (years)')
    title(['\fontsize{16} ',place(i,:),...
           '  Shift = ',num2str(shifttotal(i))])
    pause
end

%  save results 

%  results for basisw = create_bspline_basis([0,365], 1, 1)

save betaestcell1 betaestcell
save yfd1         yfd
save shift1       shifttotal
save deformmat1   deformtotal

%  results for basisw = create_bspline_basis([0,365], 5, 4)

save betaestcell5 betaestcell
save yfd5         yfd
save shift5       shifttotal
save deformmat5   deformtotal


figure(4)
subplot(1,2,1)
load betaestcell1
load yfd1
i=5;
xvec = eval_fd(daytime, tempfd(i));
yvec = eval_fd(daytime, yfd(i));
effecti = zmat(i,:);
zonei   = find(effecti == 1);
betafdPari = betaestcell{zonei};
betafdi    = getfd(betafdPari);
betaveci   = eval_fd(daytime, betafdi);
phdl=plot(daytime, xvec, 'k--', ...
          daytime, yvec, 'k-', ...
          daytime, betaveci, 'k:');
set(phdl, 'LineWidth', 1.5)
axis([0,365,-30,20])
axis('square')
xlabel('\fontsize{13} Day')
ylabel('\fontsize{13} Temperature (deg C)')
title('\fontsize{13} Shift only')
subplot(1,2,2)
load betaestcell5
load yfd5
i=5;
xvec = eval_fd(daytime, tempfd(i));
yvec = eval_fd(daytime, yfd(i));
effecti = zmat(i,:);
zonei   = find(effecti == 1);
betafdPari = betaestcell{zonei};
betafdi    = getfd(betafdPari);
betaveci   = eval_fd(daytime, betafdi);
phdl=plot(daytime, xvec, 'k--', ...
          daytime, yvec, 'k-', ...
          daytime, betaveci, 'k:');
set(phdl, 'LineWidth', 1.5)
axis([0,365,-30,20])
axis('square')
xlabel('\fontsize{13} Day')
title('\fontsize{13} Shift + Warp')



