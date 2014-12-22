% ============
% Description:
% ============
% 
% This is the main function to implement the Functional Additive Model (FAM), where the
% predictor is a function X(t_x) and the response can be either a function
% Y(t_y) or a scalar. 
% 
% Reference: M\"uller, H.G., Yao, F. (2008). Functional additive models.
% J. American Statistical Association 103, 1534?1544.
% 
% 
% IMPORTANT: This model is intended primarily for densely sampled trajectories
% 
% Includes the following steps:
% 
% 1) FPCA using the PACE method for X(t_x) and Y(t_y) (in the case that the response is functional)
% 
% 2) Computation of the additive model components
% 
% 3) Prediction of response functions for given predictor values
% 
% 4) Computation of (generalized) R-squared.
% 
% ========
% Usage:
% ========
% 
% [res, xx, yy] = FPCfam(x, t_x, y, t_y, param_X, param_Y, K_x, K_y, GCV, kern, isNewSub)
% 
% =======
% Input:
% =======
% 
% x    :   1*n cell array for predictor function x, where x{i} is the row vector of
%          measurements for the ith subject, i=1,...,n. It may contain data for subjects
%          that are used for prediction; this is controlled by "isNewSub" which is either a vector consisting of
%          0's and 1's according to whether subject is used for prediction (0) or estimation (1), or is controlled 
%          by a positive integer nn. In this case, nn is the number of subjects to be used for estimation and 
%          n-nn is the number of remaining subjects to be used for prediction, corresponding to 
%          the last n-nn data rows. When "isNewSub" is set to [], all n subjects
%          are used for estimation and no prediction will be calculated; see "isNewSub" for more details.
% 
% t_x  :   1*n cell array, t_x{i} is the row vector of time points for the ith
%          subject at which corresponding measurements x{i} are taken,
%          i=1,...,n. It contains subjects that are used for prediction.
%          See above for 2 different cases of "isNewSub" and the definition of
%          "isNewSub" for more details.
% 
% y    :   i) When no prediction is requested, that is, isNewSub = [] or
%             isNewSub = 0:
%             1*n cell array for response function y, y{i} is the response
%             row vector for the ith subject, i = 1,...,n, 
%             or a 1*n vector for scalar response y, then y(i) is the response value
%             for the ith subject, i = 1,...,n.
% 
%         ii) When prediction is requested, that is isNewSub is either a vector of
%             0's (subjects used for estimation) and 1's (subjects used for prediction but not estimation)
%             or a positive integer:
%             Cell array of dim 1*nn (functional response case) or a vector
%             of dim 1*nn (scalar response case). In either case, nn is the number of subjects 
%             used for estimation.
% 
%        iii) When prediction is requested and y is 1*n, it will be
%             truncated to 1*nn according to the isNewSub definition below.
%             (see "isNewSub" for more details).
% 
% t_y  :   1*n cell array, t_y{i} is the row vector of time points for the ith
%          subject at which corresponding measurements, y{i} are taken (for each response trajectory),
%          i = 1,...,n, or [] if y is a scalar response. See above for two different cases of 
%          "isNewSub" and the definition of "isNewSub" for more details.
% 
% param_X: an object that is returned by setOptions() that sets the input
%          arguments for FPCA() of the X (predictor) functions (for default, set param_X = []).
%          The default method for choosing the number of principal components
%          is 'AIC_R', which uses the regression AIC criterion, selecting the number
%          of principal components based on the linear relationship between predictor
%          and response (Note that if FIT = -1, the number of principal
%          components of X and Y are jointly selected by 'AIC_R'). 
%          If the main goal of the analysis is the estimation of the regression parameter function or surface, 
%          suggested settings are FIT = 0 and 'BIC1' for 'selection_k' in param_X.
%          For optimal prediction, it is better to use the default options.
%          For other default values, see setOptions().
% 
% param_Y: an object that is returned by setOptions(), that sets the input
%          arguments for FPCA() of the Y (response) functions (for default, set param_Y = []).
%          The default method for choosing the number of principal
%          components for the response trajectories is 'BIC1'. See setOptions() for more details.
%          When y is a scalar, this object will be ignored.
% 
% K_x:     positive integer, not required. Number of principal components
%          of predictor x used in regression, must be smaller than or equal to the maximum
%          number of principal components available from FPCA.  Default K_x = [], in this case
%          K_x is the number of components selected by FPCA.m for predictor functions.
% 
% K_y:     positive integer, not required. No effect for the
%          case of a scalar response. The number of principal components for the response y,
%          used for functional regression, must be smaller than or equal to the maximum
%          number of principal components available from FPCA. Default K_y = [],
%          Here K_y is the number of components selected by FPCA.m for response
%          functions.
% 
% GCV:     1  use generalized cross-validation (GCV) to choose the bandwidth when estimating the additive
%             model components (default);
%          0  use cross-validation (CV) to choose the bandwidth when estimating the additive
%             model components.
% 
% kern:    a character string for the kernel to be used when estimating the
%          additive model components
%               'epan': Epanechikov kernel
%               'gauss': Gaussian kernel  (default)
%               'gausvar': variant of Gaussian kernel 
%               'rect': rectangular kernel
%               'quar': quartic kernel
% 
% isNewSub: i) 1*n vector of 0s or 1s, where
% 
%              1 : the data for the corresponding subject, i.e.,
%                  x(isNewSub == 1), t_x(isNewSub==1), t_y(isNewSub == 1)
%                  are used for prediction only;
% 
%                 The count is n-nn for subjects with isNewSub = 1.
% 
%              0 : the data for the corresponding subject, i.e.,
%                  x(isNewSub == 0), t_x(isNewSub == 0), y and
%                  t_y(isNewSub == 0) are used for estimation only.
% 
%                  The count is nn for subjects with isNewSub = 0
% 
%           This option is convenient for computing leave-one-out
%           prediction if desired.
% 
% 
%          ii) If it is a positive integer, say 5, then the last 5 subjects (in the order from
%              top to bottom within the array x) and their values for
%              t_x, t_y are used for prediction. In other words, when choosing
%              this option, one would append the ``new'' subjects for which one desires prediction of the response
%              to occur at the end of x, t_x and t_y. Then the first nn rows of the arrays x, t_x,t_y and the entire 
%              array y (of length nn) are used for estimation and the remainder (last n-nn rows) 
%              of x,t_x and t_y (of length n-nn) will be used for prediction.
%              Example: isNewSub = 5; n = length(x)
%                       x(1:(n-5)), t_x(1:(n-5)), y, t_y(1:(n-5)) are used for estimation.
%                       x((n-4):n), t_x((n-4):n) and t_y((n-4):n) are used for prediction.
% 
%           This option is convenient for obtaining predictions for a set of new subjects.
% 
%           iii) set isNewSub = [] for the default value, which is no prediction.
% 
%        Note:
% 
%        o   nn > n - nn, the number of subjects used for estimation should generally be larger than
%            number of subjects used for prediction.
% 
%        o   When no prediction is requested, x,t_x, y, t_y will be of
%            length n = nn.
% 
%        o   When prediction is requested, x, t_x, t_y will be of length n,
%            which includes nn subjects for estimation and n-nn subjects for
%            prediction. Here y is always of length nn (< n).
% 
%        o   When prediction is requested and y is of length n < nn, a warning
%            message will be given and y will be truncated to length nn. This
%            assumes that only nn among the n available data for y are used for
%            estimation and the remaining n-nn observations will be ignored.
% 
% Details: i) There are no default values for the first 3 arguments, that
%             is, they are always part of the input arguments. Note, t_y can be
%             set to [] when y is a scalar response;
%         ii) Any unspecified or optional arguments can be set to "[]" for
%             default values;
%        iii) FPCA() calls PCA(), so setOptions() sets the input
%             arguments for PCA() and the returned object contains all
%             values returned by PCA();
%         iv) Names of objects can be retrieved through names() function i.e.,
%             names(xx), names(yy), names(res) and the actual values can be
%             retrieved through the getVal() function, example: getVal(res, 'BETA'),
%             getVal(res, 'newy') etc.
%          v) When isNewSub is set to be [], no prediction will be performed,
%             and the number of subjects in x is the same as the number of subjects
%             in y;
%             when isNewSub is either a vector of 0's and 1's or a positive
%             integer, the number of subjects in x is larger than the number of subjects
%             of y, since the former contains new subjects for prediction only,
%             whose response will not be used in the model estimation process.
%          Vi)Default: no prediction will be performed.
% 
% =======
% Output:
% =======
%   res:  an aggregated object that contains fxi, spcx, newx, new_tx, newy, new_ty, 
%         r2, fitted_y, K_x, K_y, isYFun.
% 
%         1) fxi:  a nn*K_x*K_y array containing the estimated additive
%                  model components. fxi(:,i,j) is the smoothed function
%                  between the j-th component of Y and the i-th component of X.
% 
%         2) spcx: a nn*K_y matrix containing the sorted component scores of X. 
% 
%         3) newx: 1*numNewSub cell array that contains measurements for new x
%                  (predictor) functions.
% 
%         4) new_tx: 1*numNewSub cell array that contains time points corresponding
%                    to newx.
% 
%         5) newy:   1*numNewSub cell array that contains predicted measurements for
%                    corresponding newx.
% 
%         6) new_ty: 1*numNewSub cell array that contains time points corresponding
%                    to newy.
% 
%         7) r2:     the (generalized) R-square
%                    If Y is a scalar,
%                       r2 = 1-sum((Y_i-Yhat_i)^2)/sum((Y_i-mean(Y))^2)
%                    If Y is a function,
%                       r2 = 1-\int((Y-Yhat)^2)/\int((Y-mean(Y))^2)
% 		      
%                    NOTE: This is a measure of the fraction of the variation
%                    of the responses that is explained by the regression.
%                    Similar to the Quasi R-square in FPCreg.
% 
%         9) fitted_y: 1*nn (same length of y) cell array containing fitted
%                      measurements for corresponding y values that were used for
%                      estimation; when y is scalar, it is a 1*nn vector.
%        
%         10) K_x:   number of principal components used in regression for predictor X.
%         11) K_y:   number of principal components used in regression for response Y.
% 
%  xx:   an aggregated object that contains the returned values from
%        FPCA(x,t_x,param_X).
%        See PCA() or pcaHELP.txt for more details. Use getVal() to retrieve
%        the individual values.
% 
%  yy:   an aggregated object that contains the returned values from
%        FPCA(y,t_y,param_Y) or returns "[]" if y is a scalar response
%        (see PCA() or pcaHELP.txt for more details). Use getVal() to retrieve
%        the individual values.
% 
% 
%    o    See exampleFam.m for an example of the functional response case, where
%         both predictors and responses are dense, but irregular.
% 
% See also PCA, FPCA, FPCreg

function [res, xx, yy] = FPCfam(x, t_x, y, t_y, param_X, param_Y, K_x, K_y, GCV, kern, isNewSub)


   [x, t_x, y, t_y, isYFun, newx, new_tx, newy, new_ty, invalid] = pre(x, t_x, y, t_y, isNewSub);
   nn = length(x);
   
   if invalid == 1
       return;
   end
   
   if isempty(GCV)
       GCV = 1;
   end
   if isempty(kern)
       kern = 'gauss';
   end

   if isempty(param_X)
       param_X = setOptions('selection_k','BIC1');   %set default for param_X
   elseif strcmp(getVal(param_X,'selection_k'),'BIC1') == 0
       fprintf(1,'Reminder: Suggested method of choosing the number of principal components for predictor X is BIC1!\n');
   end
   if ~isempty(K_x)
       if K_x > 0
           param_X.selection_k = K_x;
       else
           fprintf(1,'Warning: K_x must be a positive integer, reset selection_k to BIC1 now\n');
           param_X.selection_k = 'BIC1';
       end      
   end
       
   verbose_x = getVal(param_X,'verbose');
   if strcmp(verbose_x, 'on') == 1
       fprintf(1, 'Obtain functional object for x:\n');
   end
   xx = FPCA(x, t_x, param_X);   %perform PCA on x
   K_x = getVal(xx,'no_opt');
   pcx = getVal(xx,'xi_est'); 
   % the following line is added by Zhu.   
   lamhat=getVal(xx,'lambda');
   
   if isYFun == 1
       
       if isempty(param_Y)
           param_Y = setOptions('selection_k','BIC1');   %set default for param_Y
       elseif strcmp(getVal(param_Y,'selection_k'),'BIC1') == 0
           fprintf(1,'Reminder: Suggested method of choosing the number of principal components for response Y is BIC1!\n');
       end
       if ~isempty(K_y)
           if K_y > 0
               param_Y.selection_k = K_y;
           else
               fprintf(1,'Warning: K_y must be a positive integer, reset selection_k to BIC1 now\n');
               param_Y.selection_k = 'BIC1';
           end
       end

       verbose_y = getVal(param_Y,'verbose');
       if strcmp(verbose_y, 'on') == 1
           fprintf(1,'Obtain functional object for y:\n');
       end
       yy = FPCA(y, t_y, param_Y);   %perform PCA on y       
       K_y = getVal(yy,'no_opt'); 
       pcy = getVal(yy,'xi_est'); 
       muy = getVal(yy,'mu');
       phiy = getVal(yy,'phi');
       phiy = phiy(:,1:K_y);
       out1y = getVal(yy,'out1');
       
   else

      yy = [];
      pcy = y';
      K_y = 1;
      verbose_y = [];
      
   end
   
   
   if strcmp(verbose_x, 'on') == 1 || strcmp(verbose_y, 'on') == 1
       fprintf(1,'Local polynomial regression of pcy on pcx:\n');
   end
   
   bw = zeros(K_y,K_x);
   spcx = zeros(nn,K_x);
   fxi = zeros(nn,K_x,K_y);
   win = ones(1,nn);
   if GCV == 1
       for i = 1:K_x
           spcx(:,i) = sort(pcx(:,i));
           for j = 1:K_y
               bw(j,i) = gcv_lwls(pcy(:,j)',pcx(:,i)',kern,1,1,0,2,verbose_x);
               [invalid,tmp] = lwls(bw(j,i),kern,1,1,0,pcx(:,i)',pcy(:,j),win,spcx(:,i)');
               %fxi(:,i,j) = tmp';
               fxi(:,i,j) = tmp'-mean(tmp); % zhu test: adjust f_k by centering at its sample mean. Based on the assumption E(f_k)=0
               %fxi(:,i,j) = tmp'-trapz(spcx(:,i),tmp'.*normpdf(spcx(:,i),0,sqrt(lamhat(i)))); % zhu test: adjust f_k by centering at its sample mean. Based on the assumption E(f_k)=0
           end
       end       
   else
       for i = 1:K_x
           spcx(:,i) = sort(pcx(:,i));
           for j = 1:K_y
               bw(j,i) = cv_common(kern,1,1,0,pcx(:,i)',pcy(:,j));
               [invalid,tmp] = lwls(bw(j,i),kern,1,1,0,pcx(:,i)',pcy(:,j),win,spcx(:,i)');
               fxi(:,i,j) = tmp'; % note: here we should also center f_k so that E(f_k)=0.
           end
       end
   end
   
   if ~isempty(newx)
       [ypred,newpcx] = FPCApred(xx,newx,new_tx);
       clear ypred;
       newpcx = newpcx(:,1:K_x);
       nsub = size(newpcx,1);
       newpcy = zeros(nsub,K_y);
       newfxi = zeros(nsub,K_x);
       for j = 1:K_y
           for i = 1:K_x
               [~, tempidx, ~] = unique(spcx(:,i), 'first'); % This line is added by zhu. When the x in interp1(x,y,xnew) contain repetitions, error. So remove the repeted x.
               newfxi(:,i) = interp1(spcx(tempidx,i),fxi(tempidx,i,j),newpcx(:,i),'spline');
              %newfxi(:,i) = interp1(spcx(:,i),fxi(:,i,j),newpcx(:,i),'spline');
              % trainfxi(:,i) = interp1(spcx(tempidx,i),fxi(tempidx,i,j),spcx(:,i),'spline');
           end
           newpcy(:,j) = sum(newfxi,2);
           %yhat=sum(trainfxi,2);
       end
       if isYFun == 1          
           if getVal(yy,'regular') == 2
               newmuy = interp1(out1y,muy,new_ty{1},'spline');
               newphiy = interp1(out1y',phiy,new_ty{1}','spline');
               newy = mat2cell(repmat(newmuy,nsub,1)+newpcy*newphiy',ones(1,nsub))';               
           else
               newy = cell(1,nsub);
               for i = 1:nsub
                   newmuy = interp1(out1y,muy,new_ty{i},'spline');
                   newphiy = interp1(out1y',phiy,new_ty{i}','spline');
                   newy{i} = newmuy+newpcy(i,:)*newphiy';
               end
           end
       else 
           newy = mean(y)+newpcy';
       end
   end
   
   fitfxi = zeros(nn,K_x);
   fitpcy = zeros(nn,K_y);
   for j = 1:K_y
       for i = 1:K_x
           [~, tempidx, ~] = unique(spcx(:,i), 'first'); % This line is added by zhu. When the x in interp1(x,y,xnew) contain repetitions, error. So remove the repeted x.
           fitfxi(:,i) = interp1(spcx(tempidx,i),fxi(tempidx,i,j),pcx(:,i),'spline');
       end
       fitpcy(:,j) = sum(fitfxi,2);
   end   
   if isYFun == 1
       ycopy = getVal(yy,'y');
       t_ycopy = getVal(yy,'t');
       fitted_y = mat2cell(repmat(muy,nn,1)+fitpcy*phiy',ones(1,nn))';
       fitted_y1 = cell(1,nn);
       for i = 1:nn
           fitted_y1{i} = interp1(out1y,fitted_y{i},t_ycopy{i},'spline');
       end
%        if getVal(yy,'regular') == 2
%            muycopy = interp1(out1y,muy,t_ycopy{1},'spline');
%            phiycopy = interp1(out1y',phiy,t_ycopy{1}','spline');
%            fitted_y = mat2cell(repmat(muycopy,nn,1)+fitpcy*phiycopy',ones(1,nn))';
%        else
%            fitted_y = cell(1,nn);
%            for i = 1:nn
%                muycopy = interp1(out1y,muy,t_ycopy{i},'spline');
%                phiycopy = interp1(out1y',phiy,t_ycopy{i}','spline');
%                fitted_y{i} = muycopy+fitpcy(i,:)*phiycopy';
%            end           
%        end
       r2 = getQ(fitted_y1, yy, ycopy, t_ycopy, isYFun,1);
   else
       fitted_y = mean(y)+fitpcy';
       r2 = getQ(fitted_y, [], y, t_y, isYFun);
       K_y = [];
   end
        
   resNames = {'fxi', 'spcx', 'newx', 'new_tx', 'newy', 'new_ty', 'r2', 'fitted_y', 'K_x', 'K_y', 'isYFun'};
   res = {fxi, spcx, newx, new_tx, newy, new_ty, r2, fitted_y, K_x, K_y, isYFun, resNames};

end
