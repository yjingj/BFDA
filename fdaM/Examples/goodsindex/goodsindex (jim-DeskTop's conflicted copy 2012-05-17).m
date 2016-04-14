%  attach FDA functions

addpath 'c:/Program Files/Matlab/fdaM'
addpath 'c:/Program Files/Matlab/fdaM/examples/goodsindex'

%  Last modified 13 September 2009

%  ------------------------------------------------------------------------
%                    Set up the data for analysis
%  ------------------------------------------------------------------------

% %  This is input of the raw data.  See below for loading the data
%    from a .mat file.
% 
% fid = fopen('nondurprod.dat','rt');
% temp = fscanf(fid,'%f');
% temp = reshape(temp, [18, 81]);
% tempmat = temp(2:13,:);
% tempmat(12,81) = 0;
% nondurables = reshape(tempmat, [12*81, 1]);
% nondurables = nondurables(1:971);
% 
% %  for completeness, make dec 99 equal to dec 98, jan 00 equal to jan 99
% 
% nondurables = [nondurables; nondurables(961,:)];
% nondurables = [nondurables; nondurables(962,:)];
% ndur = 973;
%  
%  set up time values
% yeartime = ((0:(ndur-1))./12 + 1919)';
% 
% %  compute log nondurables
% 
% lognondur = log10(nondurables);

%  The data are in a struct object called goodsindex.  
%  This object is saved as the .mat file goodsindex.mat.  
%  The following command loads the object.

load goodsindex

ndur        = goodsindex.ndur;
yeartime    = goodsindex.yeartime;
nondurables = goodsindex.nonurables;
lognondur   = goodsindex.lognonfur;

%  compute linear trend

xmat = [ones(ndur,1), yeartime];
lognondurhat = xmat*(xmat\lognondur);

%  plot the index

plot(yeartime, nondurables, '-')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Nondurable Goods Index' )
axis([1919,2000,0,120])

%  plot the log index and its linear trend

plot(yeartime, lognondur, '-', yeartime, lognondurhat, 'g--')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Log_{10} Nondurable Goods Index' )
axis([1919,2000,0.7,2.2])

%  smooth the log data with order 8 splines, knots at data points

nbasis = 979;
norder = 8;
goodsrange = [1919,2000];
goodsbasis = create_bspline_basis(goodsrange, nbasis, norder);

Lfdobj     = int2Lfd(4);
lambda     = 10^(-11);
goodsfdPar = fdPar(goodsbasis, Lfdobj, lambda);

lognondursmth = smooth_basis(yeartime, lognondur, goodsfdPar);
              
%  plot the data and smooth for 1964-1966

index = find(yeartime >= 1964 & yeartime <= 1967);
yeartimefine = linspace(1964,1967,101)';
fit = eval_fd(yeartimefine, lognondursmth);

plot(yeartime(index), lognondur(index), 'o', ...
     yeartimefine, fit, '-')
axis([1964,1967,1.61,1.73])
hold on
plot([1965,1965],[1.61,1.73],'r:',...
     [1966,1966],[1.61,1.73],'r:')
hold off

%  estimate non-seasonal trend by smoothing with knot
%  at each year.  Use order 4 splines.

longtermbasis = create_bspline_basis(goodsrange, 83);

longtermfit = data2fd(lognondur, yeartime, longtermbasis);

%  compute and plot seasonal trend

seasonfit = lognondur - eval_fd(yeartime, longtermfit);

plot(yeartime, seasonfit, '-', goodsrange, [0,0], 'r:')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Seasonal Trend' )
axis([1919,2000,-0.08,0.08])

%  plot a sine to illustrate kinetic and potential energy

tval = (0:0.01:1)';
xval = sin(2.*pi*tval);

plot(tval, xval, '-', [0,1], [0,0], 'r:')
xlabel('\fontsize{16} Time')
ylabel('\fontsize{16} Horizontal Position')

%  a phase-plane plot of the sine

Dxval  =  (2*pi).*cos(2.*pi*tval);
D2xval = -(2*pi)^2.*sin(2.*pi*tval);

plot(Dxval, D2xval, '-', ...
     [-2*pi,2*pi], [0,0], 'r:', ...
     [0,0], [-(2*pi)^2,(2*pi)^2], 'r:')
xlabel('\fontsize{16} Velocity')
ylabel('\fontsize{16} Acceleration')
 
%  some labels for months
 
monthlabs = ['j','F','m','A','M','J','J','A','S','O','N','D'];

nyr       =  1;   %  number of startyr to plot
startyr   = 64;   %   starting year
yearindex = startyr + 1900 - 1919;  % indices of year

%  This is the code that carries out the phase-plane plotting.

index    = (1:(nyr*12+1)) + yearindex*12;
durrange = [min(yeartime(index)), max(yeartime(index))];
durfine  = linspace(durrange(1),durrange(2),401)';
durcrse  = linspace(durrange(1),durrange(2),nyr*12+1)';
D1f = eval_fd(durfine, lognondursmth, 1);
D2f = eval_fd(durfine, lognondursmth, 2);
D1c = eval_fd(durcrse, lognondursmth, 1);
D2c = eval_fd(durcrse, lognondursmth, 2);

plot([0,0], [-12,12], 'b:', [-.75,.75], [0,0], 'b:')
x0 = durrange(1);
hold on
for j = 1:nyr 
    indj = find(durfine >= x0+j-1 & durfine <= x0+j);
    if j == 1 || j > 3
        plot(D1f(indj), D2f(indj), 'b-')
    elseif j == 2
        plot(D1f(indj), D2f(indj), 'g-')
    else 
        plot(D1f(indj), D2f(indj), 'r-')
    end
    if j==1, indexc =  1:12; end
    if j==2, indexc = 13:24; end
    if j==3, indexc = 25:36; end
    if j <= 3
         for k=1:12
            text(D1c(indexc(k)), D2c(indexc(k)), monthlabs(k))
        end
    end
end
if nyr==1
    title(['\fontsize{19} Year ', num2str(durrange(1))]);
else 
    title(['\fontsize{19} Years ', num2str(durrange(1)), ...
            ' to ', num2str(durrange(2))]);
end
hold off
axis([-.75, .75, -12, 12])
xlabel('\fontsize{19} Velocity') 
ylabel('\fontsize{19} Acceleration')

%  -----------------  Figure D2D164  -----------------------

nyr       =  1;   %  number of startyr to plot
startyr   = 64;   %   starting year
yearindex = startyr + 1900 - 1919;  % indices of year

%  This is the code that carries out the phase-plane plotting.

index    = (1:(nyr*12+1)) + yearindex*12;
durrange = [min(yeartime(index)), max(yeartime(index))];
durfine  = linspace(durrange(1),durrange(2),401)';
durcrse  = linspace(durrange(1),durrange(2),nyr*12+1)';
D1f = eval_fd(durfine, lognondursmth, 1);
D2f = eval_fd(durfine, lognondursmth, 2);
D1c = eval_fd(durcrse, lognondursmth, 1);
D2c = eval_fd(durcrse, lognondursmth, 2);

plot([0,0], [-12,12], 'b:', [-.75,.75], [0,0], 'b:')
x0 = durrange(1);
hold on
indj = find(durfine >= x0 & durfine <= x0+1);
hdl  = plot(D1f(indj), D2f(indj), 'k-');
set(hdl, 'LineWidth', 2);
indexc =  1:12; 
for k=1:12
    plot(D1c(indexc(k)), D2c(indexc(k)), 'o')
    text(D1c(indexc(k))+0.05, D2c(indexc(k))+0.5, monthlabs(k))
end
hold off
axis([-.75, .75, -12, 12])
axis('square')
xlabel('\fontsize{19} Velocity') 
ylabel('\fontsize{19} Acceleration')

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/GoodsD2D164.ps'


