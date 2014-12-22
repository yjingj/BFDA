% function [K_x, K_y, xx, yy] = getNumK(xx, yy, y, t_y, param_X, param_Y, isYFun, method)
%=======
%Input:
%=======
% xx, yy:  The returned values from FPCA. See FPCA() for more details.
%      y:  the original input response.
%  t_y  :  The original input observation time for response.
% param_Y: an object that is returned by setOptions(), it sets the input
%          arguments for FPCA() of the X (predictor) functions.
% param_Y: an object that is returned by setOptions(), it sets the input
%          arguments for FPCA() of the Y (response) functions.
% isYFun:  a logical value, isYFun = 1 (Y is a functional response)
%                           isYFun = 0 (Y is a scalar response)
% method:  a string taking the following two values:
%            'YPC': fitting the functional linear regression;
%            'YCO': predicting with a functional linear regression. 
%          Refer to the input FIT of FPCreg for more details.
%=======
%Output:
%=======
%    K_x:  the best number of principal components selected for X.
%    K_y:  the best number of principal components selected for Y.
% xx, yy:  The updated values for xx, yy based on selected principal components.

function [K_x, K_y, xx, yy] = getNumK(xx, yy, y, t_y, param_X, param_Y, isYFun, method)
    
    pc_x = getVal(xx,'xi_est'); 
    n = size(pc_x,1);
    Kopt_x = getVal(xx,'no_opt');

    if isYFun == 1
        
        Kopt_y = getVal(yy,'no_opt');
        
        if strcmp(method,'YPC') == 1
            
            pc_y = getVal(yy,'xi_est');
            b = zeros(Kopt_y,Kopt_x);
            SSE = zeros(Kopt_y,Kopt_x);
            for j = 1:Kopt_y
                for i = 1:Kopt_x
                    b(j,i) = pc_y(:,j)'*pc_x(:,i)/(pc_x(:,i)'*pc_x(:,i));
                    SSE(j,i) = (pc_y(:,j)'-b(j,:)*pc_x')*(pc_y(:,j)'-b(j,:)*pc_x')';
                end
            end

            critval = zeros(Kopt_y,Kopt_x);
            ind = 1:Kopt_x;
            for j = 1:Kopt_y
                critval(j,:) = 2*ind+n*log(SSE(j,:)/n);
            end

            w_critval = mean(critval,1);
            K_x = find(w_critval==min(w_critval));
            K_y = Kopt_y;
            
        elseif strcmp(method,'YCO') == 1
            
            out1_y = getVal(yy,'out1copy');
            out21_y = getVal(yy,'out21copy');
            [lambda_y, phi_y, eigen_y, noeig_y] = getEigens(getVal(yy,'xcovcopy'),out1_y,out21_y,length(out21_y),1);
            if isempty(getVal(yy, 'sigma'))
                sigma = 0;
            else
                sigma = max(getVal(yy, 'sigmanew'),getVal(yy,'rho_opt'));
            end
            Y = cell2mat(y)'-interp1(out1_y,getVal(yy,'mucopy'),cell2mat(t_y),'spline')';
            
            phi_Y = phi_y(:,1:Kopt_y);
            lambda_Y = lambda_y(1:Kopt_y);
            [X W] = getXW(pc_x, phi_Y, out1_y, t_y, lambda_Y, sigma, 1);
            X = (cell2mat(X))';
            W = blkdiag(W{:});
            SSE = zeros(1,Kopt_x);
            for i = 1:Kopt_x
                XX = X(:,1:i*Kopt_y);
                b = pinv(XX'*W*XX)*XX'*W*Y;
                SSE(i) = (Y-XX*b)'*(Y-XX*b);
            end
            ind = 1:Kopt_x;
            critval = 2*Kopt_y*ind+length(Y)*log(SSE/length(Y));
            K_x1 = find(critval==min(critval));
            
            pc_X = pc_x(:,1:K_x1);
            [X W] = getXW(pc_X, phi_y, out1_y, t_y, lambda_y, sigma, 2);
            X = (cell2mat(X))';
            W = blkdiag(W{:});
            SSE = zeros(1,noeig_y);
            for i = 1:noeig_y
                XX = X(:,1:i*K_x1);
                b = pinv(XX'*W*XX)*XX'*W*Y;
                SSE(i) = (Y-XX*b)'*(Y-XX*b);
            end
            ind = 1:noeig_y;
            critval = 2*K_x1*ind+length(Y)*log(SSE/length(Y));
            K_y = find(critval==min(critval));
            
            phi_Y = phi_y(:,1:K_y);
            lambda_Y = lambda_y(1:K_y);
            [X W] = getXW(pc_x, phi_Y, out1_y, t_y, lambda_Y, sigma, 1);
            X = (cell2mat(X))';
            W = blkdiag(W{:});
            SSE = zeros(1,Kopt_x);
            for i = 1:Kopt_x
                XX = X(:,1:i*K_y);
                b = pinv(XX'*W*XX)*XX'*W*Y;
                SSE(i) = (Y-XX*b)'*(Y-XX*b);
            end
            ind = 1:Kopt_x;
            critval = 2*K_y*ind+length(Y)*log(SSE/length(Y));
            K_x = find(critval==min(critval));
                        
            yy = setVal(yy,'no_opt',K_y);
            if K_y > Kopt_y
                yy = setVal(yy,'lambda',lambda_Y);
                yy = setVal(yy,'phi',phi_Y);
                yy = setVal(yy,'phicopy',phi_Y);
                eigen_Y = eigen_y(:,1:K_y);
                yy = setVal(yy,'eigen',eigen_Y);
                yy = setVal(yy,'eigencopy',eigen_Y);
                xcovfit = zeros(size(eigen_Y,1),size(eigen_Y,1));
                for i = 1:K_y
                    xcovfit = xcovfit + lambda_Y(i)* eigen_Y(:,i)*eigen_Y(:,i)';
                end
                yy = setVal(yy,'xcovfit',xcovfit);
                yy = setVal(yy,'xcovfitcopy',xcovfit);
                [pc_Y, pc_Y_var]=getScores(y,t_y,getVal(yy,'mucopy'),phi_Y,lambda_Y,getVal(yy,'sigma'),K_y,getVal(param_Y,'error'),getVal(param_Y,'method'),getVal(param_Y,'shrink'),out1_y,getVal(yy,'regular'),getVal(yy,'rho_opt'));
                yy = setVal(yy,'xi_est',pc_Y);
                yy = setVal(yy,'xi_var',pc_Y_var);
            end

        end
                   
    else
        
        ind = 1:Kopt_x;
        pc_y = y-mean(y);
        b = zeros(1,Kopt_x);
        SSE = zeros(1,Kopt_x);
        for i = 1:Kopt_x
            b(i) = pc_y*pc_x(:,i)/(pc_x(:,i)'*pc_x(:,i));
            SSE(i) = (pc_y-b*pc_x')*(pc_y-b*pc_x')';
        end
        critval = 2*ind+n*log(SSE/n);
        K_x = find(critval==min(critval)); 
        K_y = [];
                
    end
    
    xx = setVal(xx,'no_opt',K_x);
    if strcmp(getVal(param_X,'verbose'), 'on') == 1
        fprintf(1,['Best number of principal components of X selected by AIC_R: ' num2str(K_x) '.\n']);
    end
    if strcmp(method,'YCO') == 1 && strcmp(getVal(param_Y,'verbose'), 'on') == 1
        fprintf(1,['Best number of principal components of Y selected by AIC_R: ' num2str(K_y) '.\n']);
    end

    