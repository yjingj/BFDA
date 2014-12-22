% ============
% Description:
% ============
% 
%   This is the main function to perform Functional Principal Component Analysis
%   for generalized logitudinal observations or generalized repeated measurements (GRM). 
%   The non-Gaussian observations Y, 
%   which are repeatedly collected for a sample of individuals over time are
%   modeled with latent Gaussian processes X.
% 
% Includes the following steps:
% 
% 1) Estimation of the mean function for Y and X
% 
% 2) Estimation of the covariance functions for Y and X
% 
% 3) Calculation of the number of selected components and the overdispersion parameter
% 
% 4) Computing the smoothed individual trajectories for X and Y
% 
% ========
% Usage:
% ========
% 
%     [yy] = FPCgrm(y,t,p)
% 
% ======
% Input: 
% ======
%      y:    1*n cell array, y{i} is the vector of non-Gaussian observations for the ith subject, i=1,...,n.
%      t:    1*n cell array, t{i} is the vector of time points for the ith subject on which
%            corresponding observations y{i} are taken, i=1,...,n.
%      p:    a struct obtained from setOptions.m sets the rest of arguments for FPCgrm.m
%            ex:
%            >> p = setOptions();
%            type 'names(p)' to see the option names included in p.
%            Most of the options are similar to those used in FPCA.m; please refer to the help file of PAC.m.
%            Special options for FPCgrm.m:
%                family: a character string naming the distribution of the repeated observations.
%                        Accepted families include: 'binomial', 'poisson'.
%                param:  parameter in the designated family. default is [] (no extra parameter needed).
% 
% =======
% Output:  
% =======  
%     yy:    a cell array that contains the following values:
%                no_opt:   integer, automatically or subjectively selected
%                          value of K, the number of selected components.
% 
%                sigma:    scalar, estimate of overdispersion parameter.
% 
%                xlambda:  1*K vector, estimated eigenvalues of latent Gaussian process X
%                          (variances of functional principal components scores).
% 
%                xphi:     N*K matrix, estimated principal component
%                          functions for latent Gaussian process X,
%                          valued at distinct input time points with ascending
%                          order of all subjects, corresponding to out1.
% 
%                xeigen:   ngrid*K matrix, estimated principal component functions 
%                          for latent Gaussian process X,
%                          valued at out21, ngrid of the pooled distinct time points
%                          with ascending order of all subjects,
%                          xphi is an interpolated version of xeigen at out1.
% 
%                xi_est:   n*K matrix, predictions for random coeffecients (PC
%                          scores) for n subjects.
% 
%                ymu:       1*N vector, estimated mean functions for observed non-Gaussian process Y,
%                           valued at distinct input time points (newdata = []), in ascending order from
%                           all subjects, corresponding to out1; when newdata is defined,
%                           corresponds to the time points from newdata, same as out1.
% 
%                xmu:       1*N vector, estimated mean functions for latent Gaussian process X,
%                           valued at distinct input time points (newdata = []), in ascending order from
%                           all subjects, corresponding to out1; when newdata is defined,
%                           corresponds to the time points from newdata, same as out1.
% 
%                bw_mu:     scalar(>0), automatically or subjectively selected
%                           bandwidth for smoothing mean curve.
% 
%                ycov:      ngrid*ngrid matrix, smoothed covariance surface for observed non-Gaussian 
%                           process Y,corresponding to out21.
% 
%                xcov:      ngrid*ngrid matrix, smoothed covariance surface for latent Gaussian process X,
%                           corresponding to out21.
% 
%                bw_cov:    1*2 vector(>0), automatically or subjectively selected
%                           bandwidths for smoothing covariance surface.
% 
%                xcovfit:   ngrid * ngrid matrix, fitted covariance surface for latent Gaussian process X, 
%                           based on truncated estimate of eigenvalues ("xlambda") and
%                           principal component functions ("xeigen"), corresponding to out21.
%                           
%                AIC:      1*maxk vector, AIC values obtained when choosing K from
%                          K=1 to K=maxk, where AIC(K) is the minimum. If AIC
%                          method is not applied, it is [].
% 
%                BIC:      1*maxk vector, BIC values obtained when choosing K from
%                          K=1 to K=maxk, where BIC(K) is the minimum. If BIC
%                          method is not applied, it is [].
% 
%                FVE:      1*ngrid vector of fraction of variance explained.
% 
%                x_pred:   1*n cell array, x_pred{i} is the vector of predictions for
%                          the latent Gaussian process X for the ith subject,
%                          evaluated at time points from the output grid vector "out1".
% 
%                y_pred:   1*n cell array, y_pred{i} is the vector of predictions for
%                          the observed non-Gaussian process Y for the ith subject,
%                          evaluated at time points from the output grid vector "out1".
% 
%             x_predOrig:  1*n cell array, x_predOrig{i} is the vector of predictions
%                          for the latent Gaussian process X for the ith subject,
%                          evaluated at the same time points as the input.
% 
%             y_predOrig:  1*n cell array, y_predOrig{i} is the vector of predictions
%                          for the observed non-Gaussian process Y for the ith subject,
%                          evaluated at the same time points as the input.
% 
%                out1:     1*N vector, distinct input time points with ascending
%                          order from all subjects if p.newdata = []; otherwise, it
%                          is the same as p.newdata.
% 
%                out21:    1*ngrid vector, a grid of time points for which the
%                          smoothed covariance surface assumes values, i.e.,
%                          ngrids from out1.
% 
%   To see an example, check with example_grm.m
function [yy] = FPCgrm(y,t,p)
       
    if strcmp(p.verbose, 'on') == 1
        fprintf(1,'Part I: Obtain smoothed mean functions\n');
    end
    [ymu,bw_mu,xmu,out1] = grm_mu_NP(y,t,p.family,p.param,p.kernel,p.bwmu,...
        p.bwmu_gcv,p.regular,p.verbose);
    
    if strcmp(p.verbose, 'on') == 1
        fprintf(1,'Part II: Obtain smoothed covariance surfaces\n');
    end
    [invalid, ycov, bw_cov, xcov, out21] = grm_cov_NP(y,t,ymu,xmu,p.family,p.param,...
               p.kernel,p.bwxcov,p.bwxcov_gcv,p.ntest1,p.ngrid1,p.ngrid,p.regular,p.error,out1,p.verbose);
     
    if strcmp(p.verbose, 'on') == 1
        fprintf(1,'Part III: Compute number of principal components and overdispersion parameter\n');
    end
    if invalid == 0
        [no_opt, FVE] = no_FVE(xcov, p.FVE_threshold);
    else
        no_opt = []; FVE = [];
        return;
    end
    no_optFVE = no_opt;
    no_optCopy = no_opt;

    pc_options = {'AIC1','BIC1','FVE','user'};
    AIC = []; BIC = [];
    if ischar(p.selection_k)
        k_id = strmatch(p.selection_k, pc_options,'exact');
        if isempty(k_id)
            fprintf(1,['Warning: Invalid method name for selection_k! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
            k_id = 3;
        end
    elseif isnumeric(p.selection_k) && p.selection_k > 0
        k_id = 4;
        no_opt = p.selection_k;
    else
        fprintf(1,['Warning: "selection_k" must be a positive integer! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
        k_id = 3;
    end
    
    if k_id == 3 || k_id == 4
        [xlambda, xphi, xeigen, no_opt] = getEigens(xcov,out1,out21,no_opt,1);
        [sigma] = overdisp(y,t,xmu,ymu,p.family,p.param,xlambda,xphi,out1);
    elseif k_id == 1 || k_id == 2
        [xlambda, xphi] = getEigens(xcov,out1,out21,p.maxk,1);
        [sigma] = overdisp(y,t,xmu,ymu,p.family,p.param,xlambda(1:no_opt),xphi(:,1:no_opt),out1);
        [no_optCopy1,cret_val] = noeig_param(p.selection_k,y,t,xmu,ymu,xlambda,xphi,sigma,p.family,p.param,out1);
        no_optCopy0 = no_optCopy;
        while no_optCopy1 ~= no_optCopy && no_optCopy1 ~= no_optCopy0
            no_optCopy0 = no_optCopy;
            no_optCopy = no_optCopy1;
            [sigma] = overdisp(y,t,xmu,ymu,p.family,p.param,xlambda(1:no_optCopy),xphi(:,1:no_optCopy),out1);
            [no_optCopy1,cret_val] = noeig_param(p.selection_k,y,t,xmu,ymu,xlambda,xphi,sigma,p.family,p.param,out1);
        end
        if no_optCopy1 == p.maxk
            fprintf(1,['Warning: ' pc_options{k_id}  ' cannot find the best No. of PC for maxk = ' num2str(maxk) '. Increase maxk to get better results.\n']);
        end
        if k_id == 1
            AIC = cret_val;
        else
            BIC = cret_val;
        end
        [xlambda, xphi, xeigen, no_opt] = getEigens(xcov,out1,out21,no_optCopy1,1);
    end
    
    if strcmp(p.verbose, 'on') == 1
        fprintf(1, ['Best number of principal components selected by ' pc_options{k_id} ': ' num2str(no_opt) '.\n']) ;
        if k_id ~= 3
            fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n']);
        else
            fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation (threshold = ' num2str(p.FVE_threshold) ').\n']);
        end
        fprintf(1,['FVE calculated from ' num2str(p.ngrid) ' possible eigenvalues: \n']);
        disp(FVE);
    end

    if strcmp(p.control, 'look')
       fve2plot = FVE(1:ceil(p.ngrid/2))*100;
       figure;
       plot(0:length(fve2plot), [0 fve2plot],'--ro', 'LineWidth',2,...
       'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
       xlabel('\bf{No. of Principal Components}');
       ylabel('\bf{FVE (%)}');
       title(['\bf{Fraction of variance explained by No. of PC (threshold = ' num2str(p.FVE_threshold) ') for latent Gaussian process X}'])
       hold on
       plot(linspace(0,no_optFVE,30),ones(30,1)*FVE(no_optFVE)*100,'b',...
       ones(30,1)*no_optFVE, linspace(0,FVE(no_optFVE)*100,30),'b');
       text(no_optFVE+0.2, FVE(no_optFVE)*100-10, {[' k = ' num2str(no_optFVE) ', FVE = ' num2str(roundoff(FVE(no_optFVE)*100,3)) '%'] ' (threshold choice)'});
       axis([0 length(FVE)+1 0 101]);
       hold off

       no_opt = input('Enter the number of principal components you want to choose:\nK=');
       if strcmp(verbose, 'on') == 1
          fprintf(1, ['You just chose ' num2str(no_opt) ' principal component(s).\n']);
          fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n\n']);
       end
    end
    
    %Now output scree plot based on final no_opt
    if p.screePlot == 1
        createSP(FVE, no_opt, 'X');
    end
    
    xcovfit = zeros(length(out21),length(out21));
    for i = 1:no_opt
        xcovfit = xcovfit+xlambda(i)*xeigen(:,i)*xeigen(:,i)';
    end
    
    if isempty(p.newdata)
        p.newdata = out1;
    end
    [xi_est,x_pred,y_pred,x_predOrig,y_predOrig] = grm_predict(y,t,xmu,ymu,xlambda,xphi,sigma,...
        [],p.family,p.param,1:length(y),out1,p.newdata);

    Ynames = {'no_opt','sigma','xlambda','xphi','xeigen','xi_est','ymu','xmu','bw_mu','ycov','xcov','bw_cov','xcovfit',...
        'AIC','BIC','FVE','x_pred','y_pred','x_predOrig','y_predOrig','out1','out21'};
    yy = {no_opt,sigma,xlambda,xphi,xeigen,xi_est,ymu,xmu,bw_mu,ycov,xcov,bw_cov,xcovfit,...
        AIC,BIC,FVE,x_pred,y_pred,x_predOrig,y_predOrig,out1,out21,Ynames};
    
    
    
    
    
    
    
    
    
    
