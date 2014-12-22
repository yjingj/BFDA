%function newy = predict(newxscore, new_tx, xx, y, yy, new_ty, b, isYFun, K_x, K_y, method, alpha)
%======
%Input:
%======
%  newxscore:  new estimated FPC scores of predictor X.
%  new_tx: 1*numNewSub cell array, new_tx{i} is the row vector of time points
%          for the ith new subject on which corresponding measurements newx{i} were taken, 
%          i = 1,...,numNewSub.
%  xx, yy: The returned values from FPCreg. See FPCreg() for more details.
%  y:      The original response variable. 
%  new_ty: 1*numNewSub cell array, new_ty{i} is the row vector of time points
%          for the ith new subject on which corresponding measurements newy{i} will be taken, 
%          i = 1,...,numNewSub, or [] if y is scalar response.
%  b:      estimated coefficients of beta in the basis represatation 
%          a K_x*K_y matrix if Y is functional; a 1*K_x vector is Y is scalar 
%  isYFun:  a logical value, isYFun = 1 (Y is a functional response)
%                           isYFun = 0 (Y is a scalar response)
%  K_x:    positive integer; number of principal components of predictor x
%          used in prediction must be smaller than or equal to that used in
%          regression.
%  K_y:    positive integer; number of principal components of response y
%          used in prediction must be smaller than or equal to that used in 
%          regression.
% method:  a string taking the following two values:
%            'YPC': fitting the functional linear regression;
%            'YCO': predicting with a functional linear regression. 
%          Refer to the input FIT of FPCreg for more details.
% alpha:   the level of the confidence bands.  alpha = 0.05 if is left empty. 
%          No confidence bands will be produced if the inputted alpha is out of (0, 1).
%=======
%Output:
%=======
%  newy:   when y is functional, 1*numNewSub cell array, newy{i} is the
%          predicted response curve for the i-th new subject;
%          when y is scalar, 1*numNewSub vector, newy(i) is the
%          predicted response for the i-th new subject.
%  newcb:  when y is functional, 2*numNewSub cell array, newcb{1,i} is the
%          lower confidence bound for the i-th new subject, while newcb{2,i} 
%          is the upper confidence bound;
%          when y is scalar, 2*numNewSub vector, newcb(1,i) is the
%          lower confidence bound for the i-th new subject, while
%          newcb(2,i) is the upper confidence bound.

function [newy newcb] = predict(newxscore, new_tx, xx, y, yy, new_ty, b, isYFun, K_x, K_y, method,  alpha)
    
    if nargin < 12
        alpha = -1;
    elseif isempty(alpha)
        alpha = 0.05;
    end
    calcb = (alpha>0 && alpha<1);
    newcb = [];
    nsub = size(newxscore,1);     
    
    xphi = getVal(xx, 'phi');
    xphi = xphi(:,1:K_x);
    xout1 = getVal(xx,'out1');
    lambda = getVal(xx,'lambda');
    lambda = lambda(1:K_x);
    sigma1 = max(getVal(xx,'sigma'),getVal(xx,'rho_opt'));
    
    if isYFun == 1
        
        b = b(1:K_x,1:K_y);
        yphi = getVal(yy, 'phi');
        yphi = yphi(:,1:K_y);
        yout1 = getVal(yy,'out1');
        ymu = getVal(yy,'mu');
        
        if strcmp(method,'YPC') == 1
                       
            if getVal(yy,'regular') == 2
                new_ty = new_ty{1};
                newymu = interp1(yout1,ymu,new_ty,'spline');
                newyphi = interp1(yout1',yphi,new_ty','spline');
                newymat = repmat(newymu,nsub,1) + newxscore * b * newyphi';
                newy = mat2cell(newymat, ones(1,nsub))';
                if calcb
                    if getVal(xx,'regular') == 2
                        tx = new_tx{1};
                        xphi = interp1(xout1',xphi,tx','spline');
                        H = repmat(lambda,length(tx),1).* xphi;
                        Omega = diag(lambda)-H'*pinv(xphi*diag(lambda)*xphi'+sigma1*eye(length(tx)))*H;
                        cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(newyphi*b'*Omega*b*newyphi'))';
                        newcbmat = [newymat - repmat(cb_width, nsub,1) newymat + repmat(cb_width, nsub,1)];
                        newcb = mat2cell(newcbmat, ones(1,nsub), length(new_ty)*ones(1,2))';
                    else
                        newcb = cell(2,nsub);
                        for j=1:nsub
                            txi = new_tx{j};
                            xphii = interp1(xout1',xphi,txi','spline');
                            H = repmat(lambda,length(new_tx),1).* xphii;
                            Omega = diag(lambda)-H'*pinv(xphii*diag(lambda)*xphii'+sigma1*eye(length(txi)))*H;
                            cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(newyphi*b'*Omega*b*newyphi'))';
                            newcb{1,j} = newy{1,j} - cb_width;
                            newcb{2,j} = newy{1,j} + cb_width;
                        end
                    end
                end
            else
                newy = cell(1,nsub);
                newcb = cell(2,nsub);
                for j = 1:nsub
                    newtyi = new_ty{j};
                    newymu = interp1(yout1,ymu, newtyi, 'spline');
                    newyphi = interp1(yout1',yphi,newtyi','spline');
                    newy{j} = newymu + newxscore(j,:) * b * newyphi';
                    if calcb
                        txi = new_tx{j};
                        xphii = interp1(xout1',xphi,txi','spline');
                        H = repmat(lambda,length(txi),1).* xphii;
                        Omega = diag(lambda)-H'*pinv(xphii*diag(lambda)*xphii'+sigma1*eye(length(txi)))*H;
                        cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(newyphi*b'*Omega*b*newyphi'))';
                        newcb{1,j} = newy{1,j} - cb_width;
                        newcb{2,j} = newy{1,j} + cb_width;
                    end
                end
            end
                        
        elseif strcmp(method,'YCO') == 1
            
            b = b';
            [X] = getXW(newxscore, yphi, yout1, new_ty, [], [], 1);
            newy = cell(1,nsub);
            for j = 1:nsub
                newy{j} = interp1(yout1,ymu,new_ty{j},'spline')+(X{j}'*b(:))';
            end
            
            if calcb
                if getVal(yy,'regular') == 2&& getVal(xx,'regular') == 2
                    tx = new_tx{1};
                    xphi = interp1(xout1',xphi,tx','spline');
                    newyphi = interp1(yout1',yphi,new_ty{1}','spline');
                    H = repmat(lambda,length(tx),1).* xphi;
                    Omega = diag(lambda)-H'*pinv(xphi*diag(lambda)*xphi'+sigma1*eye(length(tx)))*H;
                    cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(newyphi*b'*Omega*b*newyphi'))';
                    newymat = cell2mat(newy')';
                    newcbmat = [newymat - repmat(cb_width, nsub,1) newymat + repmat(cb_width, nsub,1)];
                    newcb = mat2cell(newcbmat, ones(1,nsub), length(new_ty)*ones(1,2));
                else
                    newcb = cell(2,nsub);
                    for j=1:nsub
                        txi = new_tx{j};
                        xphii = interp1(xout1',xphi,txi','spline');
                        newyphi = interp1(yout1',yphi,new_ty{j}','spline');
                        H = repmat(lambda,length(txi),1).* xphii;
                        Omega = diag(lambda)-H'*pinv(xphii*diag(lambda)*xphii'+sigma1*eye(length(txi)))*H;
                        cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(newyphi*b'*Omega*b*newyphi'))';
                        newcb{1,j} = newy{1,j} - cb_width;
                        newcb{2,j} = newy{1,j} + cb_width;
                    end
                end
            end
        end
        
    else
        
        b = b(1:K_x);
        newy = mean(y)+b*newxscore';
        
        if calcb
            if getVal(xx,'regular') == 2
                tx = new_tx{1};
                xphi = interp1(xout1',xphi,tx','spline');
                H = repmat(lambda,length(tx),1).* xphi;
                Omega = diag(lambda)-H'*pinv(xphi*diag(lambda)*xphi'+sigma1*eye(length(tx)))*H;
                cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(b*Omega*b'));
                newcb = [newy - repmat(cb_width, 1, nsub); newy + repmat(cb_width, 1, nsub)];
            else
                newcb = zeros(2,nsub);
                for j=1:nsub
                    txi = new_tx{j};
                    xphii = interp1(xout1',xphi,txi','spline');
                    H = repmat(lambda,length(txi),1).* xphii;
                    Omega = diag(lambda)-H'*pinv(xphii*diag(lambda)*xphii'+sigma1*eye(length(txi)))*H;
                    cb_width = norminv(1-alpha/2,0,1) * sqrt(diag(b*Omega*b'));
                    newcb(1,j) = newy(j) - cb_width;
                    newcb(2,j) = newy(j) + cb_width;
                end
            end
        end
        
    end
