%  Add paths to data and functions

addpath ('..\..')

%  Last modified 24 September 2009

%  ----------------------------------------------------------------------
%      First step:  input and smooth the data
%  ----------------------------------------------------------------------

%  input the data from a .csv file.  The result is a data frame

[SeabirdData, SeabirdText, SeabirdRaw] = xlsread('SeabirdData.csv');

SeabirdN = size(SeabirdData,1);

Bay = SeabirdText(2:SeabirdN+1,22);

%  select only the data for sites Uyak and Uganik, which have data 
%  from 1986 to 2005, except for 1998

index = strcmp(Bay, 'Uyak') | strcmp(Bay, 'Uganik');

SeabirdTrunc = SeabirdData(index,:);

%  select the columns with counts and the years on which counts were taken

Seabird = SeabirdTrunc(:,2:17);

%  drop birds 5 and 8 because most data for 04,05 are missing for these

birdindex = [1:4, 6:7, 9:15];
nbird     = length(birdindex);

%  extract the labels for the birds to be used

birdlabels = SeabirdText(1,2:16);
birdlabels = birdlabels(birdindex);

%  select the transect numbers

transect = unique(SeabirdTrunc(:,19));

%  indices for years except for 1998, for which all data are missing

no98 = [1:12,14:20];  
birdtime = no98 - 1;
birdRng  = [0,19];

%  set up a basis for smoothing the counts.  
%  this basis is for order 4 (cubic) B-splines, with a knot at
%  each year, except for 1998 when there were no observations

nbirdbasis = 21;
knots      = [0:11,13:19];
birdbasis  = create_bspline_basis(birdRng, nbirdbasis, 4, knots);

%  set up arrays for accumulated coefficients and couts

birdcoefmat  = zeros(nbirdbasis,2*nbird);
birdcountmat = zeros(20,2*nbird);
birdnmat     = birdcountmat;

%  set up fdPar object for smoothing

lambda = 1e1;
birdfdPar = fdPar(birdbasis, 2, lambda);

%  accumulate coefficients and counts for Uyak

for i=1:79
    %  select the data and count years for this transect
    index    = find(SeabirdTrunc(:,19) == transect(i)); 
    Seabirdi = Seabird(index,[birdindex,16]);
    timei    = Seabirdi(:,nbird+1);
    %  loop through the birds to be used, smoothing each in turn
    for j=1:nbird
        %  select times corresponding to non-missing count data
        indexj        = ~isnan(Seabirdi(:,j));
        onesj         = ones(19,1);
        nanindex      = find(isnan(Seabirdi(:,j)));
        if ~isempty(nanindex)
            onesj(nanindex) = 0;
        end
        timeij        = timei(indexj)-1986;
        Seabirdij     = Seabirdi(indexj,j);
        %  smooth the data, selecting the functional data object to output
        Seabirdfdij   = smooth_basis(timeij, Seabirdij, birdfdPar);
        birdcoefmat(:        ,j) = ...
        birdcoefmat(:        ,j) + getcoef(Seabirdfdij);
        birdcountmat(timeij+1,j) = ...
        birdcountmat(timeij+1,j) + Seabirdij; 
        birdnmat(no98,j)         = ...
        birdnmat(no98,j)         + onesj; 
    end
end

%  accumulate coefficients and counts for Uganik

for i=80:129
    %  select the data and count years for this transect
    index    = find(SeabirdTrunc(:,19) == transect(i)); 
    Seabirdi = Seabird(index,[birdindex,16]);
    timei    = Seabirdi(:,nbird+1);
    %  loop through the birds to be used, smoothing each in turn
    for j=1:nbird
       %  select times corresponding to non-missing count data
        indexj      = ~isnan(Seabirdi(:,j));
        onesj       = ones(19,1);
        onesj(isnan(Seabirdi(:,j))) = 0;
        timeij      = timei(indexj)-1986;
        Seabirdij   = Seabirdi(indexj,j);
        %  smooth the data, selecting the functional data object to output
        Seabirdfdij   = ...
            smooth_basis(timeij, Seabirdij, birdfdPar);
        birdcoefmat(:        ,j+nbird) = ...
        birdcoefmat(:        ,j+nbird) + getcoef(Seabirdfdij);
        birdcountmat(timeij+1,j+nbird) = ...
        birdcountmat(timeij+1,j+nbird) + Seabirdij;
        birdnmat(        no98,j+nbird) = ...
        birdnmat(        no98,j+nbird) + onesj; 
    end
end

%  normalize coefficients and counts by dividing by number of transects

ind1 = 1:nbird;
ind2 = ind1 + nbird;
birdcoefmat(:,ind1)  = birdcoefmat(:,ind1)/79;
birdcoefmat(:,ind2)  = birdcoefmat(:,ind2)/50;
birdcountmat(:,ind1) = birdcountmat(:,ind1)/79;
birdcountmat(:,ind2) = birdcountmat(:,ind2)/50;

%  get total counts

birdtotalcount = birdcountmat(:,ind1) +  birdcountmat(:,ind2);

plot(birdtime+1986, log10(birdtotalcount([1:12,14:20],:)), 'bo-')
xlabel('\fontsize{13} Year') 
ylabel('\fontsize{13} log10(Mean count)')

%  replace 0 counts in 1998 by NA's

birdcountmat(13,:) = nan;

%  set up the functional data object for the smoothed counts

birdfdnames{1} = 'Year'; 
birdfdnames{2} = 'Birds and sites'; 
birdfdnames{3} = 'Mean count per transect';
birdfd      = fd(birdcoefmat, birdbasis, birdfdnames);

plot(birdfd)  

%  evaluate smoothing curves at 101 points for plotting

plotyear = linspace(0,19,101)';
birdmat  = eval_fd(plotyear, birdfd);

%  plot the results for each bird in turn

for i=1:nbird
    %  set limits for ordinate
    birdtot = [birdmat(:,i);          birdmat(:,i+nbird); ...
                birdcountmat(no98,i); birdcountmat(no98,i+nbird)];
    ylimit=[0,max(birdtot)];
    %  plot the smooth curve values for Uyak
    plot(plotyear+1986, birdmat(:,i),              'b-',  ...
         plotyear+1986, birdmat(:,i+nbird),        'b--', ...
         birdtime+1986,birdcountmat(no98,i      ), 'b+', ...
         birdtime+1986,birdcountmat(no98,i+nbird), 'bo', ...
         [1986,2005], [0,0], 'r:')
    xlabel('\fontsize{13} Year')
    ylabel('\fontsize{13} Mean count per transect') 
    title(['\fontsize{13} ', birdlabels(i), ...
           ' on Uyak (-) and Uganik (--)'])
    axis([1986,2005,ylimit])
    pause
end

%  ----------------------------------------------------------------------
%    Second step:  now we prepare to do some functional data analyses
%    Because counts tend to change on a percentage basis in most 
%    situations, and also because counts tend to have standard deviations
%    roughly proportional to means, we now work with log counts
%  ----------------------------------------------------------------------

%  compute the log count data (base 10)

%  set negative values for COME bird on Uganik to 0.1

birdmat(birdmat(:,16) <= 0.1,16) = 0.1;

%  smooth the logged smooth counts

logbirdfd = smooth_basis(plotyear, log10(birdmat), birdfdPar);
logbirdfdnames{1} = 'Year'; 
logbirdfdnames{2} = 'Birds and sites'; 
logbirdfdnames{3} = 'Log mean count per transect';
logbirdfd = putnames(logbirdfd, logbirdfdnames);

%  plot the log counts for each site

plot(logbirdfd)

plot(logbirdfd(1:nbird))
title('Log counts for Uyak')

plot(logbirdfd((1:nbird)+nbird))
title('Log counts for Uganik')

%  do a principal components analysis of the log counts
%  using two principal components or harmonics

nharm  = 2;
pcastr = pca_fd(logbirdfd, nharm);

%  plot the first two harmonics:
%    the first  harmonic reflects overall level
%    the second harmonic reflects increase or decrease

plot_pca_fd(pcastr);

%  plot harmonics with a bit more refinement

birdharmfd = pcastr.harmfd;

birdharmmat = eval_fd(plotyear, birdharmfd);

meanlogbirdfd = mean(logbirdfd);
birdmeanmat   = eval_fd(plotyear, meanlogbirdfd);

subplot(1,2,1)
plot(plotyear+1986, birdmeanmat, 'b-', ...
     plotyear+1986, birdmeanmat+0.2*birdharmmat(:,1), 'b+', ...
     plotyear+1986, birdmeanmat-0.2*birdharmmat(:,1), 'b-')
     xlabel('\fontsize{13} Year') 
     ylabel('\fontsize{13} Mean count per transect')
     title('\fontsize{13} Harmonic I (86.5%)')
     axis([1986,2005,0,0.6])
     axis('square')

subplot(1,2,2)
plot(plotyear+1986, birdmeanmat, 'b-', ...
     plotyear+1986, birdmeanmat+0.2*birdharmmat(:,2), 'b+', ...
     plotyear+1986, birdmeanmat-0.2*birdharmmat(:,2), 'b-')
     xlabel('\fontsize{13} Year') 
     ylabel('\fontsize{13} Mean count per transect')
     title('\fontsize{13} Harmonic II (12.1%)')
     axis([1986,2005,0,0.6])
     axis('square')

%  these two account for 98.6% of the variance

subplot(1,1,1)
plot(1:nbirdbasis, pcastr.values, 'bo-') 
xlabel('\fontsize{13} Eigenvalue number') 
ylabel('\fontsize{13} Eigenvalue')

disp(pcastr.varprop)  %  print proportions of variance

%  plot principal component scores 

scores = pcastr.harmscr;

plot(scores(:,1), scores(:,2), 'bo', ...
     [-3.5,6], [0,0],      'b--', ...
     [0,0],    [-1.5,3.5], 'b--')
xlabel('\fontsize{13} Harmonic 1') 
ylabel('\fontsize{13} Harmonic 2')
title('\fontsize{13} Uyak + Uganik')
axis([-3.5,6,-1.5,3.5])
for i=1:length(ind1)
    text(scores(ind1(i),1)+0.1, scores(ind1(i),2), birdlabels(i))
    text(scores(ind2(i),1)+0.1, scores(ind2(i),2), birdlabels(i))
end

%  -----------------------------------------------------------------
%         Analysis of variance of the two-way table
%  -----------------------------------------------------------------

%  ---------------------------------------------------------------
%        Use linear model equations to do ANOVA
%  ---------------------------------------------------------------

%  set up design matrix using 4 dummy variables for:
%    1. Intercept, or the effect of being a fish feeder on Uyak
%    2. Crustacean/Mollusc feeder effect using fish feeding as
%       bird factor baseline
%    3. Uganik effect using Uyak as site factor baseline
%    4. Product of these effects as interaction effect
%    

Zmat = zeros(26,4);

%  Intercept or baseline effect

Intercept = ones(26,1);
Zmat(:,1) = Intercept; 

%  Crustacean/Mollusc feeding effect
 
fooddummy = zeros(13,1);
fooddummy([1,2,5,8,12,13]) = 1;
fooddummy = [fooddummy; fooddummy];
Zmat(:,2) = fooddummy;

%  Uganik effect

sitedummy = zeros(26,1);
sitedummy(14:26) = 1;
Zmat(:,3) = sitedummy;

%  Interaction effect

Zmat(:,4) = fooddummy.*sitedummy;

%  compute coefficients for regression functions

ycoefmat = getcoef(logbirdfd);
bcoefmat = Zmat\ycoefmat';
betafd   = fd(bcoefmat', birdbasis);

%  plot regression functions

plot(betafd)

betamat = eval_fd(plotyear, betafd);

%  plot intercept or constant term

plot(plotyear+1986, betamat(:,1), 'b-')
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Intercept')

%  plot effect of being a crustacean/mollusc eater

plot(plotyear+1986, betamat(:,2), 'b-')
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Crustacean/Mollusc effect')

%  plot effect of being a Uyak site

plot(plotyear+1986, betamat(:,3), 'b-')
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Uyak effect')
 
%  plot interaction effect

plot(plotyear+1986, betamat(:,4), 'b-')
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Interaction  effect')

%  compute fit functions

yhatcoefmat = (Zmat*bcoefmat)';
yhatfd      = fd(yhatcoefmat, birdbasis);
yhatmat     = eval_fd(plotyear, yhatfd);

%  compute residuals

ymat   = eval_fd(plotyear, logbirdfd);
resmat = ymat - yhatmat;

varvec    = sum(resmat.^2,2)./(26-4);
stderrvec = sqrt(varvec);

plot(plotyear+1986, stderrvec, 'b-')
xlabel('\fontsize{13} year') 
ylabel('\fontsize{13} Standard error')
axis([1986,2005,0,1])

%  ---------------------------------------------------------------
%        Use function fRegress to do ANOVA
%  Here we determined that there was no contribution from either
%  site or food*site interaction.
%  ---------------------------------------------------------------

p = 2;
xfdcell = cell(p,1);
for j=1:p 
    xfdcell{j} = Zmat(:,j);
end

%  set up the functional parameter object for the regression fns.

betabasis = create_constant_basis(birdRng);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betabasis, 2, lambda, estimate);

betacell = cell(p,1);
betacell{1} = fdPar(birdbasis);
for j=2:p 
    betacell{j} = betafdPar;
end

% try fRegress_CV for choosing the amount of smoothing

loglam = -2.5:0.5:2.5;
SSE_CV = zeros(length(loglam),1);
betafdPari = betafdPar;
for i=1:length(loglam)
       betafdPari = putlambda(betafdPari,10^loglam(i));
       betacelli = betacell;
       for j=2:p 
           betacelli{j} = betafdPari;
       end
       SSE_CVi = fRegress_CV(logbirdfd, xfdcell, betacelli);
       SSE_CV(i) = SSE_CVi;
end

phdl=plot(loglam, SSE_CV, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log smoothing parameter')
ylabel('\fontsize{13} cross validated sum of squares')

% Now do this with the optimal lambda

betafdPar = fdPar(betabasis,2,10^(0.5),estimate);

betacell = cell(p,1);
betacell{1} = fdPar(birdbasis);
for j=2:p 
    betacell{j} = betafdPar;
end

%  compute regression coefficient functions and
%  predicted functions

fRegressCell = fRegress(logbirdfd, xfdcell, betacell);

%  plot regression functions

betanames{1} = 'Intercept';
betanames{2} = 'C/M';
betanames{3} = 'Site';
betanames{4} = 'Interaction';

betaestcell = fRegressCell{4};

for j=1:p
    subplot(2,2,j)
	betaestParfdj = betaestcell{j};
    betaestfdj    = getfd(betaestParfdj);
    betaestvecj   = eval_fd(plotyear, betaestfdj);
	plot(plotyear+1986, betaestvecj, 'b-')
    xlabel('\fontsize{13} Year') 
    ylabel('\fontsize{13} Temp.')
    title(['\fontsize{13} ',betanames{j}])
end

%  plot predicted functions

yhatfdobj = fRegressCell{5};

subplot(1,1,1)
plot(yhatfdobj)

%  compute residual matrix and get covariance of residuals

yhatmat  = eval_fd(plotyear, yhatfdobj);
ymat     = eval_fd(plotyear, logbirdfd);
rmat     = ymat - yhatmat;
SigmaE   = var(rmat);

%  plot covariance surface for errors

contour(SigmaE) 
xlabel('\fontsize{13} Year') 
ylabel('\fontsize{13} Year')

surf(SigmaE)

%  plot standard deviation of errors

stddevE = sqrt(diag(SigmaE));
plot(plotyear+1986, stddevE, 'b-')
xlabel('\fontsize{13} Year') 
ylabel('\fontsize{13} Standard error')
axis([1986,2005,0,1])

%  Repeat regression, this time outputting results for
%  confidence intervals

birdbasismat = eval_basis(plotyear, birdbasis);

birdcrossprod = birdbasismat'*birdbasismat;
y2cMap = birdcrossprod\bbirdbasismat;

stderrCell = fRegress_stderr(fRegressCell, y2cMap, SigmaE);

betastderrCell = stderrCell{1};

%  plot regression function standard errors

%  plot regression functions with confidence limits

for j=1:p
    subplot(2,2,j)
	betafdParj  = betaestCell{j};
	betafdj     = getfd(betafdParj);
	betaj       = eval_fd(plotyear, betafdj);
	betastderrj = eval_fd(plotyear, betastderrCell{j});
	plot(plotyear + 1986, ...
         [betaj, betaj+2*betastderrj,  betaj-2*betastderrj], 'b', ...
         [1986,2005], [0,0], 'b:')
    xlabel('\fontsize{13} Year') 
    ylabel('\fontsize{13} Reg. Coeff.')
    title(['\fontsize{13} ',betanames{j}])
end

%  -----------------------------------------------------------------
%  Use a reduced model with only a feed effect, but add bird effects
%  In this model, birds are nested within feed groups, and their
%  effects must sum to zero within each group.
%  -----------------------------------------------------------------

%  The design matrix contains an intercept dummy variable, a 
%  feed dummy variable, and dummy variables for birds, which
%  are constrained to add to zero

Zmat = zeros(28,15);

%  Intercept or baseline effect

Intercept = ones(26,1);
Zmat(1:26,1) = Intercept;

%  Crustacean/Mollusc feeding effect
 
fooddummy = zeros(13,1);
foodindex = [1,2,5,6,12,13];
fooddummy(foodindex) = 1;
fooddummy = [fooddummy; fooddummy];
Zmat(1:26,2) = fooddummy;

birddummy = diag(ones(13,1));
Zmat( 1:13,3:15) = birddummy;
Zmat(14:26,3:15) = birddummy;

%  add the constraints of zero sum for bird effects within food groups

Zmat(27,  foodindex+2 ) = 1;
Zmat(28,-(foodindex+2)) = 1;

%  add a zero function to the logbird coefficients

logbirdcoef  = getcoef(logbirdfd);
logbirdcoef0 = [logbirdcoef, zeros(nbirdbasis,2)];
logbirdfd0 = fd(logbirdcoef0,birdbasis);

p = 15;
xfdcell = cell(p,1);
for j=1:p 
    xfdcell{j} = Zmat(:,j);
end

%  set up the functional parameter object for the regression fns.

betabasis = create_constant_basis(birdRng);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betabasis, 2, lambda, estimate);

betacell = cell(p,1);
betacell{1} = fdPar(birdbasis,2,10^(-0.5));
betacell{2} = fdPar(birdbasis,2,10^(-0.5));
for j=3:p
    betacell{j} = betafdPar;
end

%  compute regression coefficient functions and
%  predicted functions

fRegressCell = fRegress(logbirdfd0, xfdcell, betacell);


loglam     = -2:0.25:0;
SSE_CV     = zeros(length(loglam),1);
betafdPari = betafdPar;
CVobs      = 1:26;
for i=1:length(loglam)
    betafdPari = putlambda(betafdPari,10^loglam(i));
    betacelli = betacell;
    for j=1:2
        betacelli{j} = betafdPari;
    end
    SSE_CV(i) = ...
        fRegress_CV(logbirdfd0, xfdcell, betacelli,CVobs);
end

subplot(1,1,1)
plot(loglam, SSE_CV, 'bo-')
xlabel('\fontsize{13} log smoothing parameter')
ylabel('\fontsize{13} cross validated sum of squares')

%  plot regression functions

betafdnames{1} = 'Intercept';
betafdnames{2} = 'C/M';

betaestcell = fRegressCell{4};

for j=1:2
    subplot(2,1,j)
	betaestParfdj = betaestcell{j};
    betaestfdj    = getfd(betaestParfdj);
    betaestvecj   = eval_fd(plotyear, betaestfdj);
	plot(plotyear+1986, betaestvecj, 'b-')
    xlabel('\fontsize{13} Year')
    ylabel('\fontsize{13} Temp.')
    title(['\fontsize{13} ',betanames{j}])
end

%  plot predicted functions

yhatfdobj = fRegressCell{5};

subplot(1,1,1)
plot(yhatfdobj)

%  compute residual matrix and get covariance of residuals

yhatmat  = eval_fd(plotyear, yhatfdobj);
ymat     = eval_fd(plotyear, logbirdfd0);
rmat     = ymat - yhatmat;
SigmaE   = var(mat);

%  plot covariance surface for errors

contour(SigmaE) 
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Year')

surf(SigmaE)

%  plot standard deviation of errors

stddevE = sqrt(diag(SigmaE));

plot(plotyear+1986, stddevE, 'b-')
xlabel('\fontsize{13} Year') 
ylabel('\fontsize{13} Standard error')
axis([1986,2005,0,1])

%  Repeat regression, this time outputting results for
%  confidence intervals

stderrCell = fRegress_stderr(fRegressCell, y2cMap, SigmaE);

betastderrCell = stderrCell{1};

%  plot regression functions with confidence limits

for j=1:2
    subplot(2,1,j)
	betafdParj  = betaestcell{j};
	betafdj     = getfd(betafdParj);
	betaj       = eval_fd(plotyear, betafdj);
	betastderrj = eval_fd(plotyear, betastderrCell{j});
	plot(plotyear+1986, ...
         [betaj, betaj+2*betastderrj, betaj-2*betastderrj], 'b-', ...
         [1986,2005], [0,0], 'b:')
    xlabel('\fontsize{13} Year') 
    ylabel('\fontsize{13} Reg. Coeff.')
    title(['\fontsize{13} ',betanames{j}])
end



