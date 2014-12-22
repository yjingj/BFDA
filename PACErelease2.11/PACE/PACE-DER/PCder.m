%===========
%Description:      
%=========== 
%
%            This is the main program to perform Functional Principal Component
%            Analysis (FPCA) via PACE for estimation of curves and their derivatives. 
%            The principal component scores can be estimated through conditional 
%            expectation or via classical integration. For the latter, one can 
%            choose a shrinkage method for estimated scores.
%
%======
%Usage:
%======
%
% function [no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,bw_mu,xcov,bw_xcov, 
%           xcovfit, AIC,BIC,FVE,y_pred, y_predOrig, out1,out21,...
%           y,t, rcov, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, ...
%           out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops]...
%           =PCder(y,t,bwmu,bwmu_gcv, bwxcov,bwxcov_gcv, ntest1,ngrid1,selection_k, 
%                  FVE_threshold, maxk,control,regular,error, ngrid,method,shrink,...
%                  newdata, kernel, npoly, nder, method_int, numBins, yname, screePlot, designPlot, corrPlot, rho, verbose)
%
%================
%Input Arguments:
%================
%
%Input y:          1*n cell array, y{i} is the vector of measurements for the 
%                  ith subject, i=1,...,n.
%
%Input t:          1*n cell array, t{i} is the vector of time points for the 
%                  ith subject for which corresponding measurements y{i} are 
%                  available, i=1,...,n.
%
%Input bwmu:       a vector of non-negative bandwidths for mu(t) and/or its
%                  derivatives
%                  0:  use cross-validation or generalized cross-validation
%                      to choose bandwidth automatically. [Default]
%                  bwmu(j) > 0: user-specified bandwidth.
%                  * bwmu(1) is always the bandwidth for mean function with nder(1) = 0;
%                  * bwmu(j), where j > 1 is bandwidth for derivative of mean function
%                    with order nder(j).
%                  See the description of "nder" below.
%
%Input bwmu_gcv:   a vector of 0's and 1's for mu(t) and/or its derivatives
%                  For choice bwmu(j) = 0, two types of cross-validations can 
%                  be performed: One-curve-leave-out cross-validation (CV) 
%                  or  generalized cross-validation (GCV)
%                  0: CV   (may be time-consuming) 
%                  1: GCV  (faster)                [Default]                
%
%Input bwxcov:     a matrix of bandwidths for covariance surface used for 
%                  smoothing of cov(X(t), X(s)) and/or its partial derivatives
%                  w.r.t. t.
%                  each row is a 1*2 vector
%                  bwxcov(i): i-th coordinate of bandwidth vector, i=1,2.
%                  bwxcov(1)==0 & bwxcov(2)==0: use cross-validation (CV)
%                  or generalized cross-validation (GCV) for automatic
%                  selection.                      [Default]
%                  bwxcov(1)>0 & bwxcov(2)>0: user-specified bandwidths.
%                  * bwxcov(1,1) and bwxcov(1,2) are always the bandwidths choices
%                    for covariance function
%                  * bwxcov(j,1) and bwxcov(j,2), where j > 1 are bandwidths for partial derivative
%                    of the covariance function, where bwxcov(j,1) is in the derivative 
%                    direction and bwxcov(j,2) is in the smoothing direction with order nder(j). 
%                  See the description of "nder" below.
%                  
%Input bwxcov_gcv: a vector of 0's and 1's for bandwidth choice used for 
%                  smoothing of cov(X(t), X(s)) and/or its partial derivatives
%                  w.r.t. t. 
%                  If setting bwxcov = [0 0], automatic bandwidth selection
%                  by CV or GCV choices
%                  0: CV method (may be time-consuming, can be accelerated by
%                     choosing small values for ntest1 and ngrid1)..
%                  1: GCV method (faster)          [Default]
%
%Input ntest1:     integer(<=n), number of curves used for CV when choosing
%                  bandwidths for smoothing the covariance surface. The subjects
%                  in the test set are randomly selected from n subjects. Small
%                  ntest1 will accelerate CV at less accuracy. [Default is 30.]
%
%Input ngrid1:     integer, number of support points for the covariance surface
%                  in the CV procedure (selecting bandwidths of covariance).
%                  Note that input ntest1 and ngrid1 provide options to save
%                  computing time when using CV or GCV.        [Default is 30.]
%
%Input selection_k: the method of choosing the number of principal components K.
%                   'AIC1': use AIC criterion with pseudo-likelihood of
%                           measurements (marginal likelihood).
%                   'AIC2': use AIC criterion with likelihood of measurements
%                           conditional on estimated random coeffecients.
%                   'BIC1': use BIC criterion with pseudo-likelihood of
%                           measurements (marginal likelihood).
%                   'BIC2': use BIC criterion with likelihood of measurements
%                           conditional on estimated random coeffecients.
%                   'FVE' (fraction of variance explained) : use scree plot
%                           approach to select number of principal components),
%                                                     [Default]
%                    see "FVE_threshold" below.
%                    'AIC_R': use regression AIC criterion (select number
%                             of principal components based on the linear
%                             relationship between predictor and response)
%                    positive integer K: user-specified number of principal
%                                        components
%
%                   Note: BIC1 and FVE produce the most parsimonious models.
%
%Input FVE_threshold:  a positive number that is between 0 and 1 [Default is 0.85.] 
%                      It is used with the option selection_k = 'FVE' and selects 
%                      the number of principal components that explain at least 
%                      "FVE_threshold" of total variation, which is the fraction 
%                      of variance explained.
%
%Input maxk:      integer, the maximum number of principal components to consider
%                 if using automatic methods to choose K, i.e., 'AIC1', 'AIC2',
%                 'BIC1' or 'BIC2' defined by selection_k. [Default is 20.]
%                 Note, when selection_k = 'FVE' or 'AIC_R', maxk is ignored.
%
%Input control:   'auto', Select K by minimizing AIC or BIC score, or find the 
%                         first K such that the FVE_threshold is exceeded.  [Default]
%                 'look', a scree plot (FVE% Vs No. of PC) will be generated based 
%                         on K <= ngrid/2, where ngrid is the number of support points 
%                         in each direction of covariance surface. User will be prompted 
%                         to enter user-specified K after viewing scree plot. This can 
%                         can be combined with any setting of selection_k.
%
%Input regular:   0, sparse (or irregular) functional data.      
%                 1, regular data with missing values
%                 2, completely balanced (regular) data.
%                 [], automatically chosen based on the design in t.   [Default]
%
%Input error:     0, no additional measurement error assumed.       
%                 1, additional measurement error is assumed.          [Default]
%
%Input ngrid:     integer, number of support points in each direction of 
%                 covariance surface when performing principal component 
%                 analysis ( ngrid > K).
%                                                           [Default is 51.]
%Input method:    used for computing random effects \xi_{ik}                 
%                 'CE': conditional expectation method      [Default]
%                 'IN': classical integration method
%                 Note: 'CE' can be applied for sparse data or regular data, but
%                       'IN' only in the case of regular data.
%
%Input shrink:    indicator of whether applying shrinkage to estimates of random 
%                 coefficients (for regular data only)
%                 0:  no shrinkage when method = 'CE' or error = 0   [Default]
%                 1:  shrinkage when method = 'IN' and error =1, otherwise, this 
%                     will be re-set to 0.
%
%Input newdata:   a row vector of user-defined output time grids for
%                 all curves. This corresponds to "out1" in the output argument
%                 If newdata = [], then "out1" corresponds to the set of distinct
%                 time points from the pooled data.
%                 "newdata" is supposed to be a vector in ascending order on
%                  the domain of the functions.           [Default is []]
%
%Input kernel:    a 1*length(nder) cell array
%                 Each element of the cell array contains a character string to define 
%                 the kernel to be used in the 1-D or 2-D smoothing for X(t) and/or
%                 its derivatives
%                 kernel{j} = 'epan'  ==> Epanechnikov kernel [Default for dense designs with n_i >= 20]  
%                             'rect'  ==> Rectangular kernel
%                             'gauss'  ==> Gaussian kernel    [Default for sparse designs, regular designs with
%                                                              missings, dense designs for n_i < 20]
%                             'gausvar' ==> variant of Gaussian kernel
%                             'quar'  ==> Quartic kernel
%                 Note: The Gaussian kernel is overall best for sparse designs but is slower than the other kernels 
%                       and if computational speed is of concern then one may wish to use the Epanechnikov kernel 
%                       also for the case of sparse designs.
%
%Input npoly:     a vector of integers that are >= 1, when it is 1, local linear fit is used
%                 Each element in the vector, npoly(j), defines the degree of polynomial 
%                 to be used to estimate mu(t), G(t,s) and/or their derivatives
%                 Note that npoly(j) must be at least nder(j)+1 and npoly(j)- nder(j) should be
%                 an odd number, where nder is defined below. 
%                 npoly(1) is ALWAYS used for estimation of mu(t) and G(t,s). [Default: npoly = nder+1]
%
%Input nder:      a vector of integers that define the order of derivative of interest.
%                 nder(j) = 0,1,2..., etc. When nder(j) = 0, estimation of curve itself.
%                 nder(1) MUST be 0 and it is used for estimation of mu(t) and G(t,s).
%                                                                             [Default: nder = 0]
%Input method_int: integer value of 1, 2 or 3 that defines the numerical integration methods for
%                  for the estimation of derivative of eigenfunction(s)
%                  1  trapezoidal rule
%                  2  Romberg integration  [Default]
%                  3  Gaussian Quadrature
%
%Input numBins:   0: no binning
%                 a positive interger (>= 10): prebin the data with user-defined
%                 number of bins. When numBins < 10, no binning will be performed.
%		  []:  prebin the data with the following rule     [Default]
%                 
%                 i) When the input data is regular = 1 or 2
%                    m = max of n_i, where n_i is the number of repeated measurements
%                    for i-th subject.
%                 ii) regular = 0
%                    m = median of n_i
%               
%                 When m <= 20 subjects, no binning.
%                 When n <= 5000 subjects and m <= 400, no binning.
%                 When n <= 5000 subjects and m > 400, numBins = 400.
%                 When n > 5000 subjects, compute
%                 
%                 m* = max(20, (5000-n)*19/2250+400)  
%
%                 if m > m*, numBins = m*
%                 if m <= m*, no binning   
%                    
%                 This option accelerates analysis, especially for data with 
%                 many time measurements. 
%
%Input yname:     a character string which denotes the name of the current 
%                 function to be estimated.               [Default is []]
%                 It is used in the title part of the scree plot output.
%                 If it is set to [], then yname is set to be the same 
%                 name as the first input argument from PCder() or FPCder(). 
%                 For example, if this analysis is for function y,
%                 then yname = 'y'.
%
%Input screePlot: indicator of whether to create the scree plot
%                 1  a scree plot will be created     
%                 0  no scree plot will be created    [Default]    
%
%Input designPlot: indicator of whether to create the design plot
%                 1  a design plot will be created
%                 0  no design plot will be created   [Default]
%                 Interpretation of design plot: All subareas of the 
%                 domain x domain support square of the covariance surface need
%                 to be covered more or less homogeneously with pairs of design points.
%Input corrPlot:  indicator of whether to create the correlation surface plot
%                 1  a correlation surface plot will be created
%                 0  no correlation surface plot will be created    [Default]
%
%Input rho:       truncation threshold for the iterative residual that is used
%                 in the estimation of FPC scores. (see FPCscore.pdf under Help/ for more details)
%                 -1: compute unadjusted FPC scores (as in previous PACE version)
%                 >0: user-defined choice of rho
%                 0: do not set truncation threshold
%                 'cv': use randomized leave-one-measurement-out CV approach to find
%                       the optimal value of rho.     [Default]
%
%Input verbose:   a character string
%                 'on': display diagnostic messages       [Default]
%                 'off': suppress diagnostic messages
%========
% Details:
%========
%       For j = 1,...,length(nder)
%       1)       bwmu(j) = 0 and bwmu_gcv(j) = 1 ==> use GCV method to choose 
%                the bandwidth for mean function or its derivative.
%
%                For Gaussian kernel, the optimal bandwidth from GCV is multiplied
%                by a factor: a) nder(j) = 0, the factor is 1.1
%                             b) nder(j) = 1, the factor is 1.2
%                             c) nder(j) >= 2, the factor is 2.
%
%                bwmu(j) = 0 and bwmu_gcv(j) = 0 ==> use CV method to choose the 
%                bandwidth for mean function or its derivative
%
%                bwmu(j) > 0 then bwmu_gcv(j) will be ignored, subjective 
%                bandwidth choice
%
%       2)       bwxcov(j,:) = [0 0] and bwxcov_gcv(j) = 1 ==> use GCV method to 
%                choose the bandwidth for cov function or its partial derivative
%
%                For Gaussian kernel, the optimal bandwidth from GCV is multiplied
%                by a factor: a) nder(j) = 0, the factor is 1.1
%                             b) nder(j) = 1, the factor is 1.2
%                             c) nder(j) >= 2, the factor is 2.
%
%                bwxcov(j,:) = [0 0] and bwxcov_gcv(j) = 0 ==> use CV method to 
%                choose the bandwidth for cov function or its partial derivative
%
%                bwxcov(j,1) > 0 and bwxcov(j,2) > 0 then bwxcov_gcv will be 
%                ignored, subjective bandwidth choice
%
%                For nder(j) = 0, if eigenvalues estimation is of primary interest, 
%                you may want to undersmooth the covariance surface for better estimates
%                of eigenvalues.
%
%       3)       If error=1, AIC1, BIC1, AIC2, BIC2, FVE all can be used 
%                for choosing K.
%
%                If error=0, only AIC1, BIC1 and FVE can be used by 
%                definition of these criteria.
%
%       4)       If some arguments are not needed or you want to use 
%                default values, just use [] for these arguments.
%
%                Alternatively, use setOptions() to define nessary input 
%                arguments and call PCA through FPCA.
%
%       5)       In general, any output in the form of 1-dimensional 
%                vectors will be a row vector and any output related to 
%                1-dimensional cell arrays will be a row cell array.
%
%=================
%Output Arguments:
%=================
%Output no_opt:   integer, automatically or subjectively selected value of K, the number of selected components.
%
%Output sigma:    scalar, estimate of measurement error variance if 
%                 error=1, while it is [] if error=0.
%
%Output lambda:   1*K vector, estimated eigenvalues (variances of functional principal components scores).
%
%Output phi:      a 1*length(nder) cell array of estimated eigenfunctions for 
%                 X(t) and/or its derivatives.
%                 Each element of the cell array contains
%                 N*K matrix, estimated eigenfunctions 
%                 valued at distinct input time points with ascending 
%                 order of all subjects, corresponding to out1
%                 phi{1} is ALWAYS the estimated eigenfunctions
%
%Output eigen:    a 1*length(nder) cell array of estimated eigenfunctions for X(t) and/or its derivatives.
%                 Each element of the cell array contains:
%                 ngrid*K matrix, estimated eigenfunctions, 
%                 valued at out21, ngrid of the pooled distinct time points 
%                 with ascending order of all subjects,
%                 phi is an interpolated version of eigen at out1
%                 eigen{1} is ALWAYS the estimated eigenfunctions
%
%Output xi_est:   n*K matrix, predictions for random coeffecients (PC 
%                 scores) for n subjects.
%
%Output xi_var:   K*K matrix, Var(PC score)-Var(estimated PC score). The 
%                 omega matrix in equation (7) of the paper, which is used 
%                 to construct the point-wise C.I. for X_i(t)
%
%Output mu:       a 1*length(nder) cell array of estimated mu(t) and/or its derivatives.
%                 Each element of the cell array contains:
%                 1*N vector, estimated mean functions valued at distinct 
%                 input time points (newdata = []), in ascending order from 
%                 all subjects, corresponding to out1; when newdata is defined, 
%                 corresponds to the time points from newdata, same as 
%                 out1.
%                 mu{1} is ALWAYS the estimated mean function
%
%Output bw_mu:    a vector of automatically or subjectively slected
%                 bandwidth for mu(t) and/or its derivatives.
%                 Each element of the vector contains:
%                 a scalar(>0), automatically or subjectively selected 
%                 bandwidth
%
%Output xcov:     a 1*length(nder) cell array of estimated smoothed G(t,s) and/or its
%                 partial derivative w.r.t. t.
%                 Each element of the cell array contains:
%                 ngrid*ngrid matrix, smoothed covariance surface (diagnal 
%                 removed), corresponding to out21
%                 xcov{1} is ALWAYS the estimated smoothed G(t,s)
%
%Output bw_xcov:  a matrix of automatically or subjectively selected
%                 bandwidths for G(t,s). 
%                 Each row of the matrix contains:
%                 1*2 vector(>0), automatically or subjectively selected 
%                 bandwidths
%
%Output xcovfit:  a 1*length(nder) cell array of fitted covariance surface and/or its
%                 partial derivative w.r.t. t. 
%                 Each element of the cell array contains:
%                 ngrid * ngrid matrix, fitted covariance surface(or its derivatives)
%                 , based on truncated estimate of eigenvalues ("lambda") and 
%                 eigenfunctions ("eigen") or its derivatives, corresponding 
%                 to out21
%                 xcovfit{1} is ALWAYS the fitted covariance surface
%
%Output xcorr:    ngrid * ngrid matrix, fitted correlation surface, based
%                 on truncated estimate of eigenvalues ("lambda") and
%                 principal component functions ("eigen"), corresponding
%                 to out21
%
%Output AIC:      1*K vector, AIC values obtained when choosing K from 
%                 K=1 to K=maxk, where AIC(K) is the minimum. If AIC 
%                 method is not applied, it is []
%
%Output BIC:      1*K vector, BIC values obtained when choosing K from 
%                 K=1 to K=maxk, where BIC(K) is the minimum. If BIC 
%                 method is not applied, it is []
%
%Output FVE:      1*ngrid vector of fraction of variance explained
%
%Output y_pred:   a 1*length(nder) cell array of cell array for the predicted X(t) and/or
%                 its derivatives.
%                 Each element of ypred, say y_predj = y_pred{j} is a cell array, such that
%                 y_predj{i} is the vector of predictions for the ith subject evaluated at 
%                 time points from the output grid vector "out1".
%                 j = 1: curve prediction with order nder(1) = 0, 
%                 j > 1, derivative of curve prediction with order nder(j).
%
%Output y_predOrig: a 1*length(nder) cell array of cell array for the predicted X(t) and/or
%                   its derivatives.
%                   Each element of y_predOrig, say y_predOrigj = y_predOrig{j} is a cell array
%                   such that y_predOrigj{i} is the vector of predictions for the ith subject 
%                   at the same time points as the input.
%                   j = 1, curve prediction with order nder(1) = 0, 
%                   j > 1, derivative of curve prediction with order nder(j).
%
%Output out1:     1*N vector, distinct input time points with ascending 
%                 order from all subjects if newdata = []; otherwise, it 
%                 is the same as newdata.
%
%Output out21:    1*ngrid vector, a grid of time points for which the 
%                 smoothed covariance surface assumes values, i.e.,
%                 ngrids from out1.
%
%Output y:        if no binning is performed, same as the input y
%                 if binning is performed, 1 * n cell array, y{i} is a vector 
%                 of measurements after binning for subject i, i = 1,...,n
%
%Output t:        if no binning is performed, same as the input t
%                 if binning, 1 * n cell array, t{i} is a vector of
%                 time points after binning for subject i, i = 1,...,n
%                 Each value of t{i} corresponds to the midpoints of
%                 each bin.
%Output rcov:     a struct object that contains tpairn, cxxn, indx, win and cyy
%                 tpairn:  2 * N  vector denotes the  pairs of time points for subject
%                 concatenating as two vectors
%                 cxxn:    1 * N vector of raw covariance corresponding to tpairn
%                 For tpairn(1,:) == tpairn(2,:), then rcov.cxxn contains measurement
%                 errors.
%
%Output regular:  if no binning is performed or regular = 2, this is the 
%                 same as the input
%                 if binning is performed and regular = 0, it will be
%                 reset to regular = 1. In other words, after binning,
%                 the sparse and irregular case is analyzed as regular 
%                 with missings by sampled data.
%
%Output rho_opt:  if rho is set as 'cv', then rho_opt is the optimal rho obtained
%                 from the CV method, otherwise, it is the same as the input rho
%
%Output sigmanew: if rho is set as >0, 0 or 'cv', then sigmanew is the iterative
%                 residual sum of squares (see FPCscore.pdf for more details). It
%                 can be used as an estimate of the variance of the measurement
%                 errors.
%                 if rho is set as -1, then sigmanew is set to the same as output
%                 sigma. 
%
% See also FPCder, setDerOptions, showDerOptionNames, exampleDer
function [no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov,...
          xcovfit, AIC,BIC,FVE,y_pred, y_predOrig, y_predDense, out1,out21,...
	  y,t, rcov, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy,xcorr,ops]...
          =PCder(y,t,bwmu,bwmu_gcv, bwxcov,bwxcov_gcv, ntest1,ngrid1,selection_k,...
          FVE_threshold, maxk,control,regular,error, ngrid,method,shrink,...
	  newdata, kernel, npoly, nder, method_int, numBins, yname, screePlot, designPlot, corrPlot, rho,verbose)
     
     no_opt =[];sigma=[];lambda=[];phi=[];eigen=[];xi_est=[]; rcov = [];
     xi_var=[];mu=[];muDense=[];bw_mu=[];xcov=[];bw_xcov=[];xcovfit=[];AIC=[];BIC=[];FVE=[];
     y_pred =[];y_predOrig = [];y_predDense = []; out1 = [];out21 =[]; 
     rho_opt=[]; sigmanew = [];mucopy=[];phicopy=[];eigencopy=[];out1copy=[];out21copy=[];xcovcopy=[];xcovfitcopy=[];xcorr=[];ops=[];
  
cut = [];     %use default value of cut

if isempty(yname)
    yname = inputname(1);
end

if isempty(nder)
     nder = 0;
elseif nder(1) ~= 0
     fprintf(1, 'Warning: nder(1) must be 0!Reset nder = [0 nder]!\n');
     nder = [0 nder];
end

if isempty(npoly) 
     npoly = nder+1;
elseif length(npoly) ~= length(nder)
    fprintf(1, 'Warning: npoly and nder must be of the same length!Reset npoly = nder+1 now!\n');
    npoly = nder+1;
end

int_choice = [1 2 3];
if isempty(method_int)
   method_int = 2;
elseif ismember(method_int, int_choice) == 0
   fprintf(1,'method_int must be 1 (Trapezoidal rule), 2 (Romberg) or 3 (Gauss Quadrature)! Reset method_int = 2 now!\n');
   method_int = 2;
end


%Set default values for the major input arguments
%when the following input arguments are set to 
%be "[]" by the user

if isempty(bwmu)|| length(bwmu) ~= length(nder)
   bwmu = zeros(1,length(nder));    %bandwidth choice for mean function is using CV or GCV
end

if isempty(bwmu_gcv)||length(bwmu_gcv) ~= length(nder)
    bwmu_gcv = ones(1,length(bwmu)); %bandwidth choice for mean function is GCV
end

bwmuLocal = zeros(1,length(bwmu));

if isempty(bwxcov)|| size(bwxcov,1) ~= length(nder)
    bwxcov = zeros(length(bwmu),2); %bandwidth choices for covariance function is CV or GCV
end

if isempty(bwxcov_gcv)||length(bwxcov_gcv) ~= length(nder)
    bwxcov_gcv = ones(1,length(bwmu)); %bandwidth choices for covariance function is GCV
end

bwxcovLocal = zeros(1, length(bwmu));

if isempty(ntest1)
    ntest1 = 30;
end
if isempty(ngrid1)
   ngrid1 = 30;
end

if isempty(selection_k)
  selection_k = 'BIC1';
end

if strcmp(selection_k, 'FVE') && isempty(FVE_threshold)
  FVE_threshold = 0.85;
end

if isempty(error)
    error = 1;       %error assumption with measurement error
end

if isempty(maxk)
    maxk = 20;      %20 PC candidates
end

if isempty(shrink)
     shrink = 0;
end

if isempty(method)
    shrink = 0;
    method = 'CE';   %method to estimate the PC score is through conditional expectation
end

if shrink == 1 && (error ~= 1 || strcmp(method,'IN') ~= 1)
  fprintf(1,'Warning: shrinkage method only had effects when method = "IN" and error = 1 !Reset to shrink = 0 now!\n');
  shrink = 0; 
end

if isempty(control)
    control = 'auto';%automatically choose the number of PC not through interactive
                     %input by the user after visualization
end

if isempty(cut) && error == 1 
    cut = 1;         %cut off the bounary
elseif isempty(cut) 
    cut = 0;
end

if error == 0 && cut == 1
   fprintf(1,'Warning: cut is only effective if error = 1!Reset to cut = 0 now!\n')
   cut = 0;
end

if isempty(ngrid)
    ngrid = 51;      %number of output time grids is 30
end

if maxk > (ngrid-2)
    fprintf(1,'Warning: maxk can only be less than or equal to ngrid-2! Reset maxk to be ngrid now!\n');
    maxk = ngrid-2;
end

if isnumeric(selection_k)
    if selection_k > (ngrid-2)
	fprintf(1,'Warning: selection_k can only be less than or equal to ngrid-2! Reset it to be ngrid now!\n');
        selection_k = ngrid-2;
   elseif selection_k <= 0
        fprintf(1,'Warning: selection_k must be a positive integer! Reset it to BIC1 now\n');
        selection_k = 'BIC1';
        FVE_threshold = 0.85;
   end
end

if isempty(screePlot)
   screePlot = 1;
end

if isempty(designPlot)
    designPlot = 0;
end

if isempty(rho)
   rho = 'cv';
end

if isempty(verbose)
  verbose = 'on';
end

if error == 0 && (strcmp(selection_k,'AIC2') || strcmp(selection_k,'BIC2'))
     fprintf(1,'Warning: When assume no measurement error, cannot use "AIC2" or "BIC2". Reset to "BIC1" now!\n');
     selection_k = 'BIC1';
end


ncohort=length(t);     % obtain the number of curves or subjects
ni = zeros(1,ncohort);
for i = 1:ncohort
  ni(i)= length(t{i});
end

if all(ni == 1)
     fprintf(1,'Error: FPCder is aborted because the data do not contain repeated measurements!');
     return;
end



%Prebin the data if conditions are satisfied as specified in the help
%for numBins.
  
if isempty(numBins)
    [newy,newt] = binData(y,t,regular,verbose);
elseif numBins >= 10
    [newy,newt] = binData(y,t,regular,verbose,numBins);
elseif numBins == 0
    newy = [];   %no binning set by user
elseif numBins < 10 || numBins < 0
    newy = [];   %no binning due to number of bins is too small
    fprintf(1,'Warning: number of bins must be at least 10! No binning will be performed!\n');
end

if isempty(newy) == 0
  y = newy;
  t = newt;
  if regular == 0
    regular = 1;
  end
end

if isempty(regular)||~any(regular==0:2)
    regular = isregular(t);
else
    ireg = isregular(t);
    if ireg < regular
        switch ireg
            case 0
                fprintf(1,'Warning: the design is sparse but has been specified as regular or regular with missing.  No computation is performed.  Please rerun with default design setting.\n');
            case 1
                fprintf(1,'Warning: the design is regular with missing but has been specified as regular.  No computation is performed.  Please rerun with default design setting.\n');
        end
        return;
    end    
end

if isempty(kernel)
    if regular == 2 && length(t{1}) >= 20
        kernel = 'epan'; %kernel: Epanechnikov
    else
        kernel = 'gauss';   %kernel: Gaussian
    end
    kernel = cellstr(repmat(kernel,length(bwmu),1));
elseif ischar(kernel) == 1
    kernel = {kernel};
    kernel = cellstr(repmat(kernel,length(bwmu),1));
end

kernNames = {'rect','gauss','gausvar','epan','quar'};
for i = 1:length(kernel)
   if isempty(strmatch(kernel{i}, kernNames, 'exact'))
       fprintf(1,['Warning: kernel ' kernel{i} ' is unrecognizable! Reset to default now.\n']);
       if regular == 2 && length(t{1}) >= 20
           kernel = 'epan'; %kernel: Epanechnikov
       else
           kernel = 'gauss';   %kernel: Gaussian
       end
   end
end


%pool all the subjects and their corresponding time points into 1 x N vectors
%tt: vector to hold time points
%yy: vector to hold the observed measurements
%ind: indices for each subject
tt = cell2mat(t);  % 1 x N vector to hold the observed time points from all subjects
yy = cell2mat(y);  % 1 x N vector to hold the observed measurements from all subjects
%indiv= []; % 1 x N vector to hold the indices of n subjects, e.g., tt(ind == 1), yy(ind == 1) for subject 1

%for i = 1:ncohort
%     indiv = [indiv repmat(i, 1, length(t{i}))];
%end


%Initial out1 is based on the unique time points of the pooled data + the unique
%time points of "newdata", the output time grid. When newdata = [], output
%"out1" is equivalent to be the unique sorted pooled time points; otherwise, it 
%corresponds to the unique "newdata".
out1 =unique([tt newdata]);

if designPlot == 1
    createDesignPlot(t, 0, 1, 1, yname);
end

ops = struct('bwmu', bwmu, 'bwmu_gcv', bwmu_gcv, 'bwxcov', bwxcov,'bwxcov_gcv', bwxcov_gcv,...
             'ntest1', ntest1,'ngrid1', ngrid1,'selection_k',selection_k,...
             'FVE_threshold', FVE_threshold, 'maxk', maxk,'control', control,'regular', regular,...
             'error', error, 'ngrid', ngrid,'method', method,'shrink', shrink,...
             'newdata', newdata, 'kernel', kernel,'npoly', npoly, 'nder', nder, 'method_int', method_int, ...
             'numBins', numBins, 'yname', yname, 'screePlot', screePlot, 'designPlot', designPlot, 'corrPlot', corrPlot, 'rho',rho, 'verbose', verbose);

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part I: Obtain smoothed mean curve\n');
end
%when bwmu = 0 and bwmu_gcv = 0, use leave-one-curve-out CV method for bw choice
%when bwmu = 0 and bwmu_gcv = 1, use leave-one-out GCV method for bw choice
%when bwmu > 0, user-defined bw for mean function

bw_mu = bwmu;
for dg = 1:length(nder)
   if bwmu(dg)==0
      if bwmu_gcv(dg) == 1
         bw_mu(dg) = gcv_lwls(yy,tt,kernel{dg},1,npoly(dg),0,regular, verbose, bwmuLocal(dg));   %use GCV method to choose bw for mean function
         if isempty(bw_mu(dg))
	    fprintf(1,['Error: FPCder is aborted because the observed data is too sparse to estimate mean function for nder = ' num2str(nder(dg)) '!\n']);
            no_opt =[];sigma=[];lambda=[];phi=[];eigen=[];xi_est=[]; rcov = [];
            xi_var=[];mu=[];muDense=[];bw_mu=[];xcov=[];bw_xcov=[];xcovfit=[];AIC=[];BIC=[];FVE=[];
            y_pred =[];y_predOrig = [];y_predDense = [];out1 = [];out21 =[];
            rho_opt=[]; sigmanew = [];mucopy=[];phicopy=[];eigencopy=[];out1copy=[];out21copy=[];xcovcopy=[];xcovfitcopy=[];
            return;
         end
         bw_mu(dg) = adjustBW1(kernel{dg},bw_mu(dg), npoly(dg),nder(dg),regular, verbose);
      else                                                                              %use CV method to choose bw for mean function
         bw_mu(dg) = cvfda_lwls(y,t,kernel{dg},1,npoly(dg),0,regular, verbose, bwmuLocal(dg));
      end
   elseif bwmu(dg) > 0 && strcmp(verbose, 'on') == 1
      fprintf(1, ['User selected bandwidth for mu' num2str(nder(dg)) ': ' num2str(bw_mu(dg)) '\n']); 
   else
      fprintf(1,'Error: Bandwidth choice for the mean function or its derivatives must be positive!\n');
      return;
   end
end

%define the vector of case weight in the local weighted least square
%here, it is set to be one for all subjects
win1 = ones(1,length(tt));
mu = cell(1, length(bwmu));
muDense = mu;
out21 = linspace(min(out1),max(out1),ngrid);
for dg = 1:length(nder)
  [invalid, mu{dg}] = lwls(bw_mu(dg),kernel{dg},1,npoly(dg),nder(dg),tt,yy',win1,out1, bwmuLocal(dg));
  [invalid, muDense{dg}] = lwls(bw_mu(dg),kernel{dg},1,npoly(dg),nder(dg),tt,yy',win1,out21, bwmuLocal(dg));
end    

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part II: Choose bandwidth of smoothing covariance surface\n');
end

%time=cputime;
rcov = getRawCov(y,t,out1,mu{1}, regular, 0);           %obtain raw covariance;
bw_xcov = bwxcov;
for dg = 1:length(nder)
   if bwxcov(dg,1)==0 || bwxcov(dg,2)==0  
       if bwxcov_gcv == 1
         btemp = gcv_mullwlsnDer(t,ngrid1,regular,error, kernel{dg}, npoly(dg), 0, rcov, verbose, bwxcovLocal(dg)); 
         if isempty(btemp)
            fprintf(1,['Error: FPCder is aborted because the observed data is too sparse to estimate covariance function for nder = ' num2str(nder(dg)) '!\n']);
            no_opt =[];sigma=[];lambda=[];phi=[];eigen=[];xi_est=[]; rcov = [];
            xi_var=[]; xcov=[];bw_xcov=[];xcovfit=[];AIC=[];BIC=[];FVE=[];
            y_pred =[];y_predOrig = [];y_predDense=[];out1 = [];out21 =[];
            rho_opt=[]; sigmanew = [];mucopy=[];phicopy=[];eigencopy=[];out1copy=[];out21copy=[];xcovcopy=[];xcovfitcopy=[];
            return;
         else
            bw_xcov(dg,:) = btemp;
         end
         bw_xcov(dg,:) = adjustBW2(kernel{dg},bw_xcov(dg,:),npoly(dg),nder(dg),regular,verbose);
      else
         bw_xcov(dg,:) = cv_mullwlsnDer(y,t,mu{1},ntest1,ngrid1,regular,error,kernel{dg}, npoly(dg), 0 , rcov, verbose, bwxcovLocal(dg));
      end
   elseif bwxcov(dg,1) > 0 && bwxcov(dg,2) > 0
      if strcmp(verbose, 'on') == 1
          fprintf(1, ['User selected bandwidth for xcov' num2str(nder(dg)) ': (' num2str(bw_xcov(dg,1)) ',' num2str(bw_xcov(dg,2)) ')\n']);
      end
   elseif bwxcov(dg,1) < 0 || bwxcov(dg,2) < 0
      fprintf(1,'Error: Bandwidth choice for the covariance function and/or its partial derivatives must be positive!\n');
      return;
   end
end

%cputime-time
if strcmp(verbose, 'on') == 1
   fprintf(1,'Part III: Choose number of principal components functions\n');
end

AB_method = {'full','rand'};
rcov1 = rcov;
if error == 1
  tpairn = rcov1.tpairn;
  tneq=find(tpairn(1,:)~=tpairn(2,:));
  cyy = rcov1.cyy;
  rcov1.tpairn = tpairn(:,tneq);
  rcov1.cxxn=cyy(tneq);
  rcov1.win=ones(1,length(rcov1.cxxn));
  if regular == 1
     rcov1.count = rcov1.count(tneq);
  end
end

xcov = cell(1,length(bwmu));
if regular == 1
  [invalid,xcov{1}]=mullwlskDer(bw_xcov(1,:), kernel{1},rcov1.tpairn,rcov1.cxxn',rcov1.win,...
npoly(1), 0, out21,out21,bwxcovLocal(1), rcov1.count);  %smooth raw covariance;
else
  [invalid,xcov{1}]=mullwlskDer(bw_xcov(1,:), kernel{1},rcov1.tpairn,rcov1.cxxn',rcov1.win,...
 npoly(1), 0, out21,out21, bwxcovLocal(1));  %smooth raw covariance;
end

xcov{1} = (xcov{1}+xcov{1}')/2;   %transform the smoothed covariance matrix to guarantee it is a symmetric matrix.

if invalid == 0
  [no_opt, FVE] = no_FVE(xcov{1}, FVE_threshold);
else
  no_opt = []; FVE = [];
  return;
end
no_optCopy = no_opt;

pc_options = {'AIC1','AIC2', 'BIC1','BIC2','FVE','user','AIC_R'};
AIC = [];
BIC = [];

if ischar(selection_k)
    k_id = strmatch(selection_k, pc_options,'exact');
    if isempty(k_id)
       fprintf(1,['Warning: Invalid method name for selection_k! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
       k_id = 5;
    end
    if k_id == 1 || k_id == 2
      [no_opt,AIC]=no_AIC(y,t,mu{1},bw_xcov(1,:),ngrid,regular,maxk,AB_method{k_id},method, shrink, out1, out21,kernel{1}, error,cut, rcov,xcov{1},npoly(1));
    elseif k_id == 3 || k_id == 4
      [no_opt,BIC]=no_BIC(y,t,mu{1},bw_xcov(1,:),ngrid,regular,maxk,AB_method{k_id-2},method, shrink, out1, out21,kernel{1}, error,cut, rcov,xcov{1},npoly(1));
    elseif k_id == 7
       no_opt = ngrid;
    end
    if isempty(no_opt)
       k_id = 5;
       no_opt = no_optCopy;
    end
elseif isnumeric(selection_k) && selection_k > 0
    no_opt = selection_k;
    k_id = 6;
else
    fprintf(1,['Warning: "selection_k" must be a positive integer! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
    k_id = 5;
end

if k_id ~= 7
    if strcmp(verbose, 'on') == 1
       fprintf(1, ['Best number of principal components selected by ' pc_options{k_id} ': ' num2str(no_opt) '.\n']) ;
    end

    if k_id ~= 5 
         if strcmp(verbose, 'on') == 1
           fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n']);
         end
    else
         if strcmp(verbose, 'on') == 1
            fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation (threshold = ' num2str(FVE_threshold) ').\n']);
         end
    end 

    if no_opt == maxk && k_id < 5
       fprintf(1,['Warning: ' pc_options{k_id}  ' cannot find the best No. of PC for maxk = ' num2str(maxk) '. Increase maxk to get better results.\n']);
    end
  
    if strcmp(verbose, 'on') == 1
       fprintf(1,['FVE calculated from ' num2str(ngrid) ' possible eigenvalues: \n']);
       disp(FVE);
    end

    if strcmp(control, 'look')
       fve2plot = FVE(1:ceil(ngrid/2))*100;
       figure;
       plot(0:length(fve2plot), [0 fve2plot],'--ro', 'LineWidth',2,...
       'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
       xlabel('\bf{No. of Principal Components}');
       ylabel('\bf{FVE (%)}');
       title(['\bf{Fraction of variance explained by No. of PC (threshold = ' num2str(FVE_threshold) ') for function ' yname '}'])
       hold on
       plot(linspace(0,no_optCopy,30),ones(30,1)*FVE(no_optCopy)*100,'b',...
       ones(30,1)*no_optCopy, linspace(0,FVE(no_optCopy)*100,30),'b');
       text(no_optCopy+0.2, FVE(no_optCopy)*100-10, {[' k = ' num2str(no_optCopy) ', FVE = ' num2str(roundoff(FVE(no_optCopy)*100,3)) '%'] ' (threshold choice)'});
       axis([0 length(FVE)+1 0 101]);
       hold off

       no_opt = input('Enter the number of principal components you want to choose:\nK=');
       if strcmp(verbose, 'on') == 1
          fprintf(1, ['You just chose ' num2str(no_opt) ' principal component(s).\n']);
          fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n\n']);
       end
    end

    %Now output scree plot based on final no_opt
    if screePlot == 1
      createSP(FVE, no_opt, yname);
    end

end

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part IV: Perform principal components analysis\n'); 
end

%Initialize the cell arrays
phi = cell(1,length(bwmu));
eigen = phi;
xcovfit = phi;
y_pred = phi;
y_predOrig = phi;
y_predDense = phi;

if error==1
    [invalid, sigma]=pc_covE(t,out1,bw_xcov(1,:),ngrid,cut,kernel{1},rcov, npoly(1));
    if invalid==0
      [xi_est, xi_var, lambda, phi{1}, eigen{1}, no_opt, xcovfit{1},y_predOrig{1}, rho_opt,sigmanew]=pc_est(y, t, mu{1}, xcov{1}, sigma, no_opt,...
													    error, method, shrink, out1, out21,regular,rho, verbose);
        if length(nder) > 1
             [invalid, xcov(2:end), phi(2:end), eigen(2:end), xcovfit(2:end)]=pc_covEder(y,t,[],bw_xcov,lambda,eigen{1},...
ngrid,out1,kernel,rcov, bwxcovLocal, npoly, nder, method_int);
             for dg = 2:length(nder)
                y_predOrig{dg} = FPCderEval(out1, mu, phi, xi_est, no_opt, 1:ncohort,t,dg);
             end
        end
    else 
        xi_est=[]; xi_var=[]; lambda=[]; phi=[]; eigen =[]; xcovfit = []; y_predOrig = [];
        return;
    end
elseif error==0
    sigma=[];
    if invalid==0
      [xi_est, xi_var, lambda, phi{1}, eigen{1}, no_opt, xcovfit{1},y_predOrig{1}, rho_opt,sigmanew]=pc_est(y, t, mu{1}, xcov{1}, sigma, no_opt, ...
													    error, method, shrink, out1, out21,regular,rho,verbose);
        if length(nder) > 1
            [invalid, xcov(2:end), phi(2:end), eigen(2:end), xcovfit(2:end)]=pc_covNEder(y,t,[],bw_xcov,lambda,eigen{1},...
      ngrid,out1,kernel,rcov, bwxcovLocal, npoly, nder, method_int);
            for dg = 2:length(nder)
                y_predOrig{dg} = FPCderEval(out1, mu, phi, xi_est, no_opt, 1:ncohort,t,dg);
            end
        end
    else 
        xi_est=[]; xi_var=[]; lambda=[]; phi=[]; eigen = []; xcovfit = []; y_predOrig = [];
        return; 
    end
end

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part V: compute the smoothed individual trajectories\n');
end

%Save these copies for internal use in the FPCreg(), where the time points are guarantteed
%to be evaluated at the distinct time points of pooled t's.
mucopy = mu;
phicopy = phi;
eigencopy = eigen;
out21copy = out21;
out1copy = out1;
xcovcopy = xcov;
xcovfitcopy = xcovfit;

%Map the results to the user-defined output time grid "newdata"
if length(newdata) > 0
   
   out21 = linspace(min(newdata),max(newdata),ngrid);
   out1 = newdata;
   XI = repmat(out21,length(out21),1);
   XI = XI(:)';
   YI = repmat(out21,1,length(out21));

   for dg = 1:length(nder)
      mu{dg} = interp1(out1copy,mu{dg}, newdata, 'spline');
      phi{dg} = interp1(out1copy, phi{dg}, newdata, 'spline');  %eigen function based on newdata
      eigen{dg} = interp1(out21copy, eigen{dg}, out21,'spline');%eigen function based on ngrid's of newdata
      if no_opt == 1
          if size(phi{dg},2) > 1
             phi{dg} = phi{dg}';
          end
          if size(eigen{dg},2) > 1
             eigen{dg} = eigen{dg}';
          end
      end

      %covariance functions based on newdata
      xcov{dg} = reshape(interp2(out21copy,out21copy,xcov{dg},XI,YI,'spline'),length(out21),length(out21));
      xcovfit{dg} = reshape(interp2(out21copy,out21copy,xcovfit{dg},XI,YI,'spline'),length(out21),length(out21));
   end

end
 
 
if invalid==0
    for dg = 1:length(nder)
       tmp = repmat(mu{dg}, ncohort,1)+xi_est*phi{dg}';
       tmp = num2cell(tmp,2);
       y_pred{dg} = tmp';
       tmp = repmat(muDense{dg}, ncohort,1)+xi_est*eigen{dg}';
       tmp = num2cell(tmp,2);
       y_predDense{dg} = tmp';
    end
end

if error==0
    sigma=[];
end

xcorr = zeros(length(out21), length(out21));
xcovfit1 = xcovfit{1};
for i = 1:length(out21)
    for j = 1:length(out21)
        xcorr(i,j) = xcovfit1(i,j)/sqrt(xcovfit1(i,i)*xcovfit1(j,j));
    end
end
%xcorr = max(min(xcorr,1),-1);

if corrPlot == 1
    if no_opt == 1
        fprintf(1,'Warning: Correlation surface is not available when only one principal component is used.');
    else
        figure;
        mesh(out21,out21,xcorr);
        xlabel('\bf{t}');
        ylabel('\bf{t}');
        zlabel('\bf{Correlation}');
        title(['\bf{Fitted correlation surface for function ' yname '}']);
    end
end
end

function [regular] = isregular(t)
    tt = cell2mat(t);
    f = length(tt)/length(unique(tt))/length(t);
    if f==1
        regular = 2;
    elseif f>.75
        regular = 1;
    else
        regular = 0;
    end
end


