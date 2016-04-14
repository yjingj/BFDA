%  Analyses of Swedish mortality data, described in the book in Section 10.4, 
%  pages 163-165.

%  "The data are the number of deaths at each age for women born in each year
%  from 1751 to 1914 and for ages 0 to 80."  (Quoted from page 163)

%  The data actually analyzed are log hazard rates, defined as
%  " ... the natural logarithm of th ratio of the number of females who die
%    at a specific age to the number of females alive at that age." 
%  (quoted from book, page 163)

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

addpath(fdaMPath)

%  Input the log hazard rates from 1751 to 1884.

fid        = fopen('SwedeLogHazard.txt','rt');
tempvec    = fscanf(fid,'%f');
SwedeLogHazard = reshape(tempvec,144,81)';

%  Input the log hazard rates for 1914

fid        = fopen('SwedeLogHazard1914.txt','rt');
SwedeLogHazard1914 = fscanf(fid,'%f');

load SwedeLogHazard

plot(0:80,SwedeLogHazard(:,1))
line(0:80,SwedeLogHazard(:,144))

%  indices of years 1751, 1810 and 1860

threeyears = [1751, 1810, 1860] - 1750;

%  plot Figure 10.10

tempLogHazard = [SwedeLogHazard(:,threeyears), SwedeLogHazard1914];

SwedeTime = (0:80)';

phdl = plot(SwedeTime, tempLogHazard, '-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} age')
ylabel('\fontsize{13} log hazard')
legend('\fontsize{12} 1751', '\fontsize{12} 1810', ...
       '\fontsize{12} 1860', '\fontsize{12} 1914', ...
       'location', 'north')

%  plot the data for 1751 to 1884

plot(SwedeLogHazard)

save SwedeLogHazard SwedeLogHazard
save SwedeLogHazard1914 SwedeLogHazard1914

load SwedeLogHazard
load SwedeLogHazard1914

%  set up a saturated order six basis for smoothing the data

SwedeRng = [0,80];

nbasis = 85;
norder =  6;
SwedeBasis = create_bspline_basis(SwedeRng,nbasis,norder);

%  set the smoothing parameters

lambda = 1e-7;
SwedefdPar = fdPar(SwedeBasis, 4, lambda);

%  smooth the data

Swedefd  = smooth_basis(SwedeTime,SwedeLogHazard,SwedefdPar);

plotfit_fd(SwedeLogHazard, SwedeTime, Swedefd)

% set up coefficient functions

nbasis = 23;
SwedeBetaBasis = create_bspline_basis(SwedeRng, nbasis);

SwedeBeta0fdPar = fdPar(SwedeBetaBasis, 2, 1e-5);

SwedeBeta1fd    = bifd(zeros(23), SwedeBetaBasis, SwedeBetaBasis);
SwedeBeta1fdPar = bifdPar(SwedeBeta1fd, 2, 1e3);

p = 2;
SwedeBetaCell = cell(p, 1);
SwedeBetaCell{1} = SwedeBeta0fdPar;
SwedeBetaCell{2} = SwedeBeta1fdPar;

%  Define the dependent and independent variables

NextYear = Swedefd(2:144);
LastYear = Swedefd(1:143);

Swede_linmodStr = linmod(NextYear, LastYear, SwedeBetaCell);

Swede_beta1fd = eval_bifd(SwedeTime, SwedeTime, ...
                          Swede_linmodStr.beta);

% Figure 10.11

surfc(SwedeTime, SwedeTime, Swede_beta1fd)
xlabel('\fontsize{13} age')
ylabel('\fontsize{13} age')
zlabel('\fontsize{13} \beta(s,t)')



