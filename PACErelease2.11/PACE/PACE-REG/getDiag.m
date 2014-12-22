% function [h d] = getDiag(res, xx, yy, param, K_x, K_y, isYFun)
% getDiag    returns hat matrix and Cook's distance vector, and draws
%            residual diagnostic plots, leverage plots and Cook's distance
%            plot, by utilizing the output of FPCreg.
% Input:
%      res:  the output from FPCreg
%   xx, yy:  The returned values from FPCA. See FPCA() for more details.
%            Let yy=y if y is scalar.
%    param:  an object that is returned by setOptions(), it sets the input 
%            arguments for FPCA() of the fitted Y processes and the residual processes. 
%            The default method of choosing the number of principal
%            components is 'FVE'. For other default values, see setOptions() for more details. 
%            When y is a scalar, this will be ignored.
%      K_x:  number of principal components of predictor X used in regression
%      K_y:  number of principal components of response Y used in regression
% Output:
%        h: a n*n matrix; estimated hat matrix 
%        d: a 1*n matrix; estimated Cook's distance

function [h d] = getDiag(res, xx, yy, param, K_x, K_y)
    
    h=[]; d=[];
    isYFun = getVal(res, 'isYFun');
    method = getVal(res, 'method');
    if strcmp(method,'YPC') == 0 && isYFun == 1
        fprintf(1,'Warning: Diagnostics of the fitted functional linear regression model is available only for the case FIT = 0.\n');
        return;
    end
    if isempty(K_x) 
        K_x = getVal(res,'K_x');
    end
    if isempty(K_y)&&(isYFun==1) 
        K_y = getVal(res,'K_y');
    end
    
    pc_x = getVal(xx,'xi_est');
    pc_x = pc_x(:,1:K_x);
    fitted_y = getVal(res, 'fitted_y1');
    n = length(pc_x(:,1));
    h = zeros(n,n);
    d = zeros(1,n);
    if isYFun == 1
        pc_y = getVal(yy,'xi_est');
        pc_y = pc_y(:,1:K_y);
        y = getVal(yy,'y');
        t_y = getVal(yy, 't');
        lambda_y = getVal(yy,'lambda');
        lambda_y = lambda_y(1:K_y);
        b = getB(pc_x, yy, y, t_y, isYFun, K_x, K_y, method);

        resd = cell(1,n);
        for i = 1:n
            resd{i}=y{i}-fitted_y{i};
        end
        if isempty(param)
            param = setOptions('screePlot',0);
        end
        pc_r = getVal(FPCA(resd, t_y, param),'xi_est');
        pc_f = getVal(FPCA(fitted_y, t_y, param),'xi_est');
    else
        pc_y = (yy-mean(yy))';
        lambda_y = 1;
        K_y = 1;
        b = getB(pc_x, [], yy, [], isYFun, K_x, K_y, method);
        b = b';
        pc_r = (yy - fitted_y)';
        pc_f = fitted_y';
    end

    mse = zeros(K_x, K_y, n);
    for i = 1:n
        for k=1:K_y
            for j=1:K_x
                mse(j,k,i) = pc_y(i,k)-pc_x(i,j)*b(j,k);
            end
        end
    end
    mse = sum(mse.^2,3)/(n-1);

    for j = 1:K_x
        h = h + pc_x(:,j)*(pc_x(:,j)'*pc_x(:,j))^(-1)*pc_x(:,j)';
    end

    if isYFun == 1
        for i = 1:n
            yyi = setVal(yy, 'xi_est', pc_y(1:n~=i,:));
            bi = getB(pc_x(1:n~=i,:), yyi, y(1:n~=i), t_y(1:n~=i), isYFun, K_x, K_y, method);
            for k = 1:K_y
                for j = 1:K_x
                    d(i) = d(i) + lambda_y(k) * (pc_x(:,j)'*pc_x(:,j)) * (b(j,k)-bi(j,k))^2 / mse(j,k);
                end
            end
        end
    else
        for i = 1:n
            bi = getB(pc_x(1:n~=i,:), [], pc_y(1:n~=i,:)', [], isYFun, K_x, K_y, method);
            bi = bi';
            for j = 1:K_x
                d(i) = d(i) + lambda_y(k) * (pc_x(:,j)'*pc_x(:,j)) * (b(j,k)-bi(j,k))^2 / mse(j,k);
            end
        end
    end
    d = d / sum(lambda_y);

    figure;
    ylim=1.1*max(max(abs(pc_r)));
    if isYFun == 1
        mr=size(pc_r,2);
        mf=size(pc_f,2);
        for j = 1:mr
            for k = 1:mf
                subplot(mr, mf,  (mr-j)*mf+k);
                plot(pc_f(:,k), pc_r(:,j),'.');
                hold on;
                plot([min(pc_f(:,k))-0.1*range(pc_f(:,k)) max(pc_f(:,k))+0.1*range(pc_f(:,k))], [0 0], 'g-');
                axis([min(pc_f(:,k))-0.1*range(pc_f(:,k)) max(pc_f(:,k))+0.1*range(pc_f(:,k)) -ylim ylim]);
                xlabel(['No.' num2str(k) ' Fitted FPC Score']);
                ylabel(['No.' num2str(j) ' Residual FPC Score']);
            end
        end
    else
                plot(pc_f, pc_r,'.');
                hold on;
                plot([min(pc_f(:,k))-0.1*range(pc_f(:,k)) max(pc_f(:,k))+0.1*range(pc_f(:,k))], [0 0], 'g-');
                axis([min(pc_f(:,k))-0.1*range(pc_f(:,k)) max(pc_f(:,k))+0.1*range(pc_f(:,k)) -ylim ylim]);
                xlabel('Fitted');
                ylabel('Residual');
    end
    
    figure;
    subplot(1, 2, 1);
    stem(diag(h),'MarkerSize',3);
    hold on;
    plot([0 n+1],[2*K_x/n 2*K_x/n], 'g--');
    xlabel('Index');
    ylabel('Functional Leverages');
    axis([0 n+1 0 max(diag(h))*1.1]);
    subplot(1, 2, 2);
    stem(d,'MarkerSize',3);
    xlabel('Index');
    ylabel('Functional Cooks Distance');
    axis([0 n+1 0 max(d)*1.1]);
    set(gca,'YGrid','on');
end
