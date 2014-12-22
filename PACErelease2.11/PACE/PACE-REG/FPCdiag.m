% ============
% Description:
% ============
% 
% This is the driver function to perform Functional Linear Regression (FLM),
% which includes functional regression diagnostics and a boostrap test for the overall
% significance of the functional regression, where the predictor is a
% function X(t_x) and the response can
% be either a function Y(t_y) or a scalar.
% 
% This function includes the following steps:
% 
% 1) FPCA using the PACE method for X(t_x) and/or Y(t_y).
% 
% 2) Computation of the regression parameter function \beta.
% 
% 3) Prediction of response functions, given observations of predictor trajectories. 
% 
% 4) Computation of functional R-square and Quasi R-square (a simple predictive
%    measure based on the concept of variance explained by the regression).
% 
% 5) Functional regression diagnostics (described in Chiou, J.M., M\"uller,
%    H.G. (2007). Diagnostics for functional regression via residual
%    processes. Computational Statistics and Data Analysis 51, 4849-4863).
% 
% 6) Boostrap test for regression effect, for the null hypothesis
%    that the functional R-square equals zero (testing for overall functional
%    regression effect, as described in M\"uller, H.G., Chiou, J.M., Leng, X. (2008).
%    Inferring gene expression dynamics via functional regression analysis. 
%    BMC Bioinformatics 9:60).
% 
% ========
% Usage:
% ========
% 
% [res xx yy h d rej pv cb grid_x grid_y] = FPCdiag(x, t_x, y, t_y,...
%         param_X, param_Y, FIT, K_x, K_y, isNewSub, bootstrap, nsim, alpha)
% 
% =======
% Input:
% =======
% 
% x    :  1*n cell array for predictor function x, where x{i} is the row vector of
%         measurements for the ith subject, i=1,...,n. It contains those subjects
%         that are used for prediction, when "isNewSub" is defined as either
%         a vector of 0's and 1's or as a positive integer. Then, nn is the number
%         of subjects to be used for estimation and n-nn is the number of subjects
%         to be used for prediction. When "isNewSub" is set to [], all n subjects
%         are used for estimation and no prediction will be calculated.
% 
%         See "isNewSub".
% 
% t_x  :  1*n cell array, t_x{i} is the row vector of time points for the ith
%         subject at which the corresponding measurements x{i} are taken,
%         i=1,...,n. It contains subjects that are used for prediction.
%         See above for two different cases of "isNewSub" and the definition of
%         "isNewSub" for more details.
% 
% y    :  i) When no prediction is requested, that is, isNewSub = [] or
%            isNewSub = 0:
%            1*n cell array for the case of a response function y, where y{i} is the response
%            row vector for the ith subject, i = 1,...,n
%            or a 1*n vector in the case of a scalar response y, where y(i) is the response value
%            for the ith subject, i = 1,...,n.
% 
%        ii) When prediction is requested, that is isNewSub is either a vector of
%            0's and 1's or a positive integer, then it is a cell array of dim 1*nn or a vector
%            of dim 1*nn, where nn is the number of subjects used for estimation.
% 
%       iii) When prediction is requested and y is of dim 1*n, it will be
%            truncated to dim 1*nn, according to the isNewSub definition;
%            see "isNewSub" for more details.
% 
% t_y  :  1*n cell array, t_y{i} is the row vector of time points for the ith
%         subject on which corresponding measurements y{i} are taken,
%         i = 1,...,n, or [] if y is a scalar response.
%         See above for the two different options for "isNewSub", and the definition of
%         "isNewSub" for more details.
% 
% param_X: an object that is returned by setOptions() and defines the input
%         arguments for the FPCA() of the X (predictor) functions (for default, set param_X = []).
%         The default method of choosing the number of principal components
%         is 'FVE'. For the other defaults, see setOptions().
% 
% param_Y: an object that is returned by setOptions(), and defines the input
%          arguments for the FPCA() of the Y (response) functions (for default, set param_Y = []).
%          The default method of choosing the number of principal
%          components for the response trajectories is 'BIC1'. See setOptions() for more details.
%          When y is a scalar, this object will be ignored.
% 
% K_x:   positive integer, normally not needed. Number of principal components
%        of predictor x used for functional regression, must be smaller than or equal to the maximum
%        number of principal components available from FPCA. Default K_x = [], then
%        K_x is the number of components selected by FPCA.m for predictor functions.
% 
% K_y:   positive integer, normally not needed. This has no effect in the
%        the case of a scalar response. For functional response, it is the 
%        number of principal components for the response y
%        used for functional regression, and it must be smaller than or equal to the maximum
%        number of principal components available from FPCA. Default K_y = [].
%        
% 
%  FIT:    an indicator with values 0 or -1:
%             FIT =  0: Fitting a functional linear regression (default) through decomposition into
%                       simple linear regressions through the origin between the FPC scores of Y (responses)
%                       and those of X (predictors) (denoted as 'YPC' in the following).
%             FIT = -1: Prediction with a functional linear regression
%                       with a functional response Y, specifically adapted to the case of sparsely sampled response 
%                       trajectories and using a linear regression of observed responses Y on FPC scores of X, 
%                       applying weighted least squares with the covariance surface of Y as the
%                       weight matrix (denoted as 'YCO' in the following).
%             Note:   The option FIT = -1 aims at optimal
%                     prediction for sparsely observed functional responses (i.e. when 
%                     Y is sparse and functional), and does not target the regression parameter surface                     \beta. 
% 
% 
% isNewSub: i) 1*n vector of 0s or 1s, with elements:
% 
%              1   if the data for the corresponding subject, i.e.,
%                  x(isNewSub == 1), t_x(isNewSub==1), t_y(isNewSub == 1)
%                  are used for prediction;
% 
%                 There will be n-nn subjects with isNewSub = 1.
% 
%              0   if the data for the corresponding subject, i.e.,
%                  x(isNewSub == 0), t_x(isNewSub == 0), y and
%                  t_y(isNewSub == 0) are used for estimation.
% 
%                  There will be nn subjects with isNewSub = 0
% 
%              Note: This option is convenient to compute leave-one-out
%              prediction if desired.
% 
% 
%          ii) If this is a positive integer, say 5, then the last 5 subjects from the arrays x,
%              t_x, t_y are used for prediction. In other words, when choosing
%              this option, append the new subjects to the end of x, t_x and
%              t_y. Then the subjects with data in the first nn rows of arrays x, t_x,t_y and all y (of length nn)
%              are used for estimation and the subjects with data in the remaining rows of x,t_x and t_y (of length
%              n-nn) are used for prediction.
%              Example:  isNewSub = 5; n = length(x)
%                  x(1:(n-5)), t_x(1:(n-5)), y, t_y(1:(n-5)) are used for
%                  estimation, while x((n-4):n), t_x((n-4):n) and t_y((n-4):n) are used for
%                  prediction
%              This option is convenient to perform prediction for a large number of new subjects.
% 
%           iii) set isNewSub = [] for the default value, i.e., no prediction.
% 
%        Note:
% 
%        o   nn > n - nn, and should not be too small. In other words,
%            the number of subjects used for estimation should be larger than
%            the number of subjects used for prediction.
% 
%        o   When no prediction is requested, x, t_x, y, t_y all will be of
%            length of n = nn.
% 
%        o   When prediction is requested, x, t_x, t_y will be of length of n,
%            which includes nn subjects for estimation and n-nn subjects for
%            prediction. Here y is always of length nn (< n).
% 
%        o   When prediction is requested and y is of length n (> nn), a warning
%            message will be issued and y will be truncated to length nn. This
%            assumes that only nn among the n available data for y are used
%            for estimation and the remaining n-nn data will be ignored.
% 
% bootstrap: Bootstrap method.  If bootstrap = 'pre', confidence bands for beta and
%           R2test are based on bootstrapping by resampling from the original subject
%           data, which is more accurate but computationally intensive and thus very slow. The 
%           default [] bootstrap is based on resampling from the estimated functional 
%           principal scores. The variability coming from both measurement errors and errors 
%           in estimating mean and eigenfunctions is neglected in this variant, which is 
%           therefore less accurate, but much faster.
% 
%   nsim:   size of the bootstrap sample 
%              bootstrap = 'pre'  ==>  default is nsim = 200
%              bootstrap = 'post' ==>  default is nsim = 500 
% 
%   alpha:  level of hypothesis tests and confidence bands, optional; pointwise
%           100(1-alpha)% confidence bands will be created for the regression
%           parameter function/surface. The null hypothesis of no regression
%           relation is tested at level alpha, based on bootstrapping the
%           functonal R-square under the null hypothesis. The
%           default value for alpha is 0.05.
% 
% Details: i) There are no default values for the first 3 arguments, that
%             is, these data must be provided as part of the input arguments. 
%             Note t_y can be set to [] when y is a scalar response.
%         ii) Any unused or optional arguments can be set to "[]" to request
%             default values;
%        iii) FPCA() calls PCA(); setOptions() sets the input
%             arguments for PCA() and the returned object contains all
%             values returned by PCA();
%         iv) Names of objects can be retrieved through the names() function, i.e.,
%             names(xx), names(yy), names(res) and their actual values can be
%             retrieved through the getVal() function, example: getVal(res, 'BETA'),
%             getVal(res, 'newy') etc.
%         v) When isNewSub is set to be [], no prediction will be performed,
%            and the number of subjects in x is the same as the number of subjects
%            in y;
%            when isNewSub is either a vector of 0's and 1's or a positive
%            integer, then the number of subjects in x is larger than the number of subjects
%            of y, since the former contains new subjects for which no response data are used 
%            and prediction of responses is requested. Any such data will not be used in the 
%            model estimation process. By default, no prediction will be performed.
% 
% =======
% Output:
% =======
%   res:  an aggregated object that contains BETA, b, newx, newtx, newy,
%         new_ty, r2, cov_xy, sigma_xy.
% 
%         1) BETA:  an object that contains beta, grid_x and grid_y.
% 
%                   i) beta: a ngrid_x*ngrid_y matrix; the estimated regression parameter
%                            surface or function, evaluated at given grid points.
%                  ii) grid_x : 1*ngrid_x vector, contains ngrid distinct time
%                               points for the x function (predictor function)
%                               where ngrid is defined in setOptions().
%                 iii) grid_y : 1*ngrid_y vector, contains ngrid distinct time
%                               points for the y function (response function),
%                               where ngrid is defined in setOptions()
%                               when y is a scalar, it is [].
%         2)  b:  a K_x*K_y matrix of estimated coefficients for the regression parameter
%                 function in the eigenbasis representation, where K_x is
%                 the number of FPC for x and K_y is the number of FPC for y; alternatively
%                 interpreted as regression slope coefficients of the response FPC scores
%                 on predictor FPC scores.
%         3) newx: 1*numNewSub cell array contains measurements for new x
%                  (predictor) functions.
% 
%         4)new_tx: 1*numNewSub cell array containing time points corresponding
%                    to newx.
% 
%         5) newy:   1*numNewSub cell array containing predicted measurements for
%                    corresponding newx.
% 
%         6) new_ty: 1*numNewSub cell array containing time points corresponding
%                    to newy.
% 
%         7) r2:     functional R-square.
% 
%         8) Q:      Quasi R-square; a measure of the fraction of variance of the
%                    responses that is explained by the regression:
% 
%                    If Y is a scalar,
%                       Q = 1-sum((Y_i-Yhat_i)^2)/sum((Y_i-mean(Y))^2).
%                    If Y is a function,
%                       Q=1-sum((Y_i-Yhat_i)'*(Y_i-Yhat_i)/n_i)/
%                            sum((Y_i-mean(Y))'*(Y_i-mean(Y))/n_i).
% 
%  xx:   an aggregated object that contains the returned values from
%        FPCA(x,t_x,param_X).
%        See PCA() or pcaHELP.txt for more details. Use getVal() to retrieve
%        the individual values.
% 
%  yy:   an aggregated object that contains the returned values from
%        FPCA(y,t_y,param_Y) or returns [] if y is a scalar response
%        (see PCA() or pcaHELP.txt for more details). Use getVal() to retrieve
%        the individual values.
% 
%  h:    a n*n matrix; estimated functional hat matrix.
% 
%  d:    a 1*n matrix; estimated functional Cook's distances.
% 
%  rej:  1=reject the null hypothesis of no regression relation; 0 = do not reject
%        (depending on alpha).
% 
%  pv:   p-value for the null hypothesis of no regression relation.
% 
%  cb:   a ngrid_x*ngrid_y*2 matrix;
%        cb(:,:,1) is the estimated lower confidence band at given grid points;
%        cb(:,:,2) is the estimated upper confidence band at given grid points;
% 
%  grid_x: 1*ngrid_x vector, contains ngrid distinct time points for the x (predictor)
%          function, where ngrid is defined in setOptions().
% 
%  grid_y: 1*ngrid_y vector, contains ngrid distinct time points for the y (response)
%          function, where ngrid is defined in setOptions(), or returns [] if y is a scalar.
% 
%  See
%    o    example0.m for an example of sparse irregular data,
%         where both predictors and responses are functional
%    o    example2.m for an example of regular data,
%         where both predictors and responses are functional
%    o    example_scal for an example of regular data,
%         where the predictor is functional and the response is scalar.
% 
% See also FPCA, FPCreg, getDiag, getCB and Rtest
% Diagnostics of the fitted functional linear regression model is available
% only for the case FIT = 0.


function [res xx yy h d rej pv cb grid_x grid_y Qsim] = FPCdiag(x, t_x, y, t_y, param_X, param_Y, FIT, K_x, K_y, isNewSub, bootstrap, nsim, alpha)
    
    if nargin < 15 || isempty(alpha)
        alpha = 0.05;
    end
    
   if isempty(FIT) || FIT == 0
       FIT = 0;
   elseif checkY(y) == 1
       fprintf(1,'Warning: Diagnostics of the fitted functional linear regression model is available only for the case FIT = 0.\n');
       fprintf(1,'Now running FPCdiag with option FIT = 0.\n');
       FIT = 0;
   end
    
    if strcmp(bootstrap,'pre')==1
        [res, xx, yy] = FPCreg(x, t_x, y, t_y,  param_X, param_Y, FIT, K_x, K_y, isNewSub, alpha);
        BETA = getVal(res, 'BETA');
        grid_x = getVal(BETA, 'grid_x');
        grid_y = getVal(BETA, 'grid_y');
        beta = getVal(BETA, 'beta');
        [x, t_x, y, t_y, isYFun] = pre(x, t_x, y, t_y, isNewSub);
                
        if isYFun == 0
            yy=y;
        end
                
        if isempty(param_X)
            param_X = setOptions();
        end
        if isempty(param_Y)
            param_Y = setOptions();
        end
        param_X.screePlot = 0;
        param_Y.screePlot = 0;

        [h d] = getDiag(res, xx, yy, param_Y, [], []);
        n=length(y);
        
        if isempty(nsim)
            nsim = 200;
        end
        rnd = mysample(1:n, n*nsim*3,1);
        rnd = reshape(rnd, n, nsim, 3);
        idx = ones(nsim, length(grid_x));
        
        figure;
        if isYFun == 1
            beta2=zeros(length(grid_x),length(grid_y),nsim);
            idy = ones(nsim, length(grid_y));
            for l=1:nsim
                rnd1 = rnd(:, l, 1);
                [res1] = FPCreg(x(rnd1), t_x(rnd1), y(rnd1), t_y(rnd1),  param_X, param_Y, FIT, K_x, K_y, 0);
                BETA1 = getVal(res1, 'BETA');
                beta1 = getVal(BETA1, 'beta');
                grid_x1 = getVal(BETA1, 'grid_x');
                idx(l,:) = (grid_x<=grid_x1(end) & grid_x>=grid_x1(1));
                grid_y1 = getVal(BETA1, 'grid_y');
                idy(l,:) = (grid_y<=grid_y1(end) & grid_y>=grid_y1(1));
                beta1 = interp2(grid_x1,grid_y1,beta1,grid_x(idx(l,:)==1),grid_y(idy(l,:)==1)','spline');
                beta2(idx(l,:)==1,idy(l,:)==1,l) = beta1;
            end
            bselx=find(sum(idx,1)==nsim);
            bsely=find(sum(idy,1)==nsim);
            beta2=beta2(bselx,bsely,:);
            S=length(bselx);
            T=length(bsely);
            cb=zeros(S,T,2);
            for j=1:S
                for k=1:T
                    cb(j,k,:)=quantile(beta2(j,k,:),[alpha/200 1-alpha/200]);
                end
            end
            
        surf(grid_x(bselx),grid_y(bsely),cb(:,:,1));
        hold on;
        surf(grid_x(bselx),grid_y(bsely),cb(:,:,2));
        surf(grid_x(bselx),grid_y(bsely),beta(bselx,bsely));
        hold off;
        xlabel('x');
        ylabel('y');
        title('Confidence bands for Beta');
        else
            beta2=zeros(length(grid_x),nsim);
            for l=1:nsim
                rnd1 = rnd(:, l, 1);
                [res1] = FPCreg(x(rnd1), t_x(rnd1), y(rnd1), [],  param_X, param_Y, FIT, K_x, K_y, 0);
                BETA1 = getVal(res1, 'BETA');
                beta1 = getVal(BETA1, 'beta');
                grid_x1 = getVal(BETA1, 'grid_x');
                idx(l,:) = (grid_x<=grid_x1(end) & grid_x>=grid_x1(1));

                beta1 = interp1(grid_x1,beta1,grid_x(idx(l,:)==1),'spline');
                beta2(idx(l,:)==1,l) = beta1;
            end
            bselx=find(sum(idx,1)==nsim);
            S=length(bselx);
            cb =zeros(S,2);
            for j=1:S
                cb(j,:)=quantile(beta2(j,:),[alpha/2 1-alpha/2]);
            end

            plot(grid_x, cb(:,1),'--');
            hold on;
            plot(grid_x, cb(:,2),'--');
            plot(grid_x, beta(bselx));
            hold off;
            xlabel('x');
            title('Confidence bands for Beta');
        end
        
        Qsim=zeros(nsim,1);

        for l=1:nsim
            rndx = rnd(:, l, 2);
            rndy = rnd(:, l, 3);
            if isYFun==1
                [res1] = FPCreg(x(rndx), t_x(rndx), y(rndy), t_y(rndy),  param_X, param_Y, FIT, K_x, K_y, 0);
            else
                [res1] = FPCreg(x(rndx), t_x(rndx), y(rndy), [],  param_X, param_Y, FIT, K_x, K_y, 0);
            end
            Qsim(l) = getVal(res1, 'Q');
        end

        pv = sum(Qsim>getVal(res,'Q'))/nsim;
        if pv < alpha
            rej = 1;
            if pv==0
                fprintf(1,['P-value < ' num2str(1/nsim) '. Reject the null hypothesis of no regression relation.\nX and Y are significantly correlated.\n']);
            else
                fprintf(1,['P-value = ' num2str(pv) '. Reject the null hypothesis of no regression relation.\nX and Y are significantly correlated.\n']);
            end
        else
            rej = 0;
            fprintf(1,['P-value = ' num2str(pv) '. Cannot reject the null hypothesis of no regression relation.\nX and Y are not significantly correlated.\n']);
        end
    else
        isYFun = checkY(y);
        [res, xx, yy] = FPCreg(x, t_x, y, t_y,  param_X, param_Y, FIT, K_x, K_y, isNewSub, alpha);
        
        if isYFun == 0
            yy=y;
        end

        [h d] = getDiag(res, xx, yy, param_Y, [], []);
        [cb grid_x grid_y] = getCB(res, xx, yy, [], [], alpha);
        [rej pv Qsim] = Rtest(res, xx, yy,nsim, alpha);
    end
end
