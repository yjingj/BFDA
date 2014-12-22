%============
%Description:
%============
%
% This is the main function to perform Functional Singular Component
% Analysis of two functions X(t_x) and Y(t_y).
%
% Includes the following steps:
%
% 1) Apply the functional singular value decomposition to the kernel
%    estimator of cross-covariance surface between X and Y.
%
% 2) Computation of the singular values and the singular functions.
%
% 3) Computation of the functional correlation.
%
%========
% Usage:
%========
%
% res = FSVD(x, t_x, param_X, y, t_y, param_Y, nsvd, bwccov)
%
%=======
% Input:
%=======
%
% x    :  1*n cell array for function x,  x{i} is the row vector of measurements
%         for the ith subject, i=1,...,n. 
%
% t_x  :  1*n cell array, t_x{i} is the row vector of time points for the ith
%         subject at which corresponding measurements x{i} are taken,
%         i=1,...,n.
%
% y    :  1*n cell array for function y,  y{i} is the row vector of measurements
%         for the ith subject, i=1,...,n. 
% 
% t_y  :  1*n cell array, t_y{i} is the row vector of time points for the ith
%         subject at which corresponding measurements y{i} are taken,
%         i=1,...,n.
%
% param_X: an object that is returned by setOptions(), it sets the input
%         arguments for FPCA() of the X (predictor) functions (for default, set param_X = []).
%         For default values of parameters, see setOptions() for more details.
%
% param_Y: an object that is returned by setOptions(), it sets the input
%         arguments for FPCA() of the X (predictor) functions (for default, set param_X = []).
%         For default values of parameters, see setOptions() for more details.
%
% nsvd:   positive integer. Number of singular components used in regression.
%         Default nsvd = [], then use 'FVE' (fraction of variance explained) 
%         criteria to select the number of singular components that explain at least
%         75% of total variation.
%
% bwccov:  1*2 vector, bandwidths for covariance surface used for
%          smoothing of cov(X(t),Y(s))
%          bwccov(i): ith coordinate of bandwidth vector, i=1,2.
%          bwccov(1)>0 & bwccov(2)>0: user-specified bandwidths.
%          bwccov(1)==0 & bwccov(2)==0: use generalized cross-validation
%          (GCV) for automatic selection.                     [Default]
%          For the purpose of estimating singular value, the GCV-chosen 
%          bandwidth is multiplied by an empirical factor depending on the
%          sparsity of the data.
%          
%
% Details: i) Any non-used or optional arguments can be set to "[]" for
%             default values;
%         ii) FPCA() calls PCA(), so setOptions() sets the input
%             arguments for PCA() and the returned object contains all
%             values returned by PCA();
%        iii) Names of objects can be retrieved through names() function i.e.,
%             names(xx), names(yy), names(res) and the actual values can be
%             retrieved through getVal() function, example: getVal(res, 'BETA'),
%             getVal(res, 'newy') etc.
%
%=======
% Output:
%=======
%   res:  an aggregated object that contains nsvd, fc, sv, ss_x, ss_y, sc_x, sc_y, out_x, out_y, out2x, out2y, ccov_s, ccovfit.
%
%         1) nsvd: automatically or subjectively selected value of the number of singular components.
%
%         2) fc:   functional correlation, a scalar between 0 and 1 .
%
%         3) sv:   a 1*nsvd vector, estimated singular values (covariances of
%                  functional singluar components scores).
%
%         4) ss_x: n*nsvd matrix, predictions for random coeffecients
%                  (singular scores) of X for n subjects.
%
%         5) ss_y: n*nsvd matrix, predictions for random coeffecients
%                  (singular scores) of Y for n subjects.
%
%         6) out_x: 1*Nx vector, distinct input time points of X with ascending
%                   order from all subjects.
%
%         7) out_y: 1*Ny vector, distinct input time points of Y with ascending
%                   order from all subjects.
% 
%         8) sc_x:  Nx*nsvd matrix, estimated singular component functions
%                   valued at distinct input time points with ascending
%                   order of all subjects, corresponding to out_x.
%
%         9) sc_y:  Ny*nsvd matrix, estimated singular component functions
%                   valued at distinct input time points with ascending
%                   order of all subjects, corresponding to out_y.
%
%        10) out2x, out2y: both are 1*51 vectors, provides a grid of time points 
%                   for which the smoothed cross-covariance surface assumes values.
%
%        11) ccov_s: 51*51 matrix, smoothed covariance surface, corresponding to 
%                    out2x and out2y.
%
%        12) ccovfit:  51*51 matrix, fitted cross-covariance surface, based
%                   on truncated estimate of singular values and component 
%                   functions, corresponding to out2x and out2y.
%
%
%  See
%    o    example0.m for an example of sparse irregular data case,
%         when both functions are functional
%    o    example2.m for example of regular data case,
%         when both functions are functional
%
% See also PCA, FPCA

function [res invalid]=FSVD(x,t_x,p_x,y,t_y,p_y,nsvd,bwccov)

if isempty(p_x)||isempty(p_y)
   regular = 0;
else
    regular = min(getVal(p_x,'regular'),getVal(p_y,'regular'));
end
if nargin<8
    bwccov = [];
end
error = ones(1,2);
if isempty(p_x)||getVal(p_x,'error')==1
    error(1)=1;
else
    error(1)=0;
end
if isempty(p_y)||getVal(p_y,'error')==1
    error(2)=1;
else
    error(2)=0;
end
% n = length(x);

[res1 invalid] = getSC(x,t_x,p_x, y,t_y,p_y, bwccov,nsvd,regular);
sv = getVal(res1,'lambda');
lambda_x = getVal(res1,'lambda_x');
lambda_y = getVal(res1,'lambda_y');
nsvd = getVal(res1, 'nsvd');

if invalid==0
    [ss_x, ss_y, ss_var, rho_opt, sig1]=getSScores(res1, x, t_x, y, t_y, nsvd, error, 'CE', 1,  regular, 'cv');
end
fc = sv(1)/sqrt(lambda_x(1)*lambda_y(1));

resnames = ['fc' 'sv' 'ss_x' 'ss_y' 'rho_opt' 'sig1' getVal(res1,'names')];
res = [fc sv ss_x ss_y rho_opt sig1 res1];
res{end} = resnames;
fprintf(1, 'The functional correlation between X and Y is: %f\n', fc);
end
