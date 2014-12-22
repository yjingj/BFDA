function [X W] = getXW(pc_x, phi_y, out1y, t_y, lambda_y, sigma, method)
% method: different ways of calculating X (but results are the same)

    [n,Kx] = size(pc_x);
    Ky = size(phi_y,2);
    X = cell(1,n);
    
    if isempty(lambda_y) && isempty(sigma)       % unweighted method
        
        W = [];
        if method == 1
            for i = 1:n
                phitmp = interp1(out1y',phi_y,t_y{i}','spline');
                pctmp = repmat(pc_x(i,:),Ky,1);
                pctmp = diag(pctmp(:));
                X{i} = pctmp*repmat(phitmp,1,Kx)';
            end
        elseif method == 2
            for i = 1:n
                phitmp = interp1(out1y',phi_y,t_y{i}','spline');
                phitmp = reshape(repmat(phitmp,Kx,1),length(t_y{i}),Kx*Ky);
                pctmp = repmat(pc_x(i,:),Ky,1)';
                pctmp = diag(pctmp(:));
                X{i} = pctmp*phitmp';
            end
        end
        
    else
        
        W = cell(1,n);
        if method == 1
            for i = 1:n
                phitmp = interp1(out1y',phi_y,t_y{i}','spline');
                W{i} = pinv(phitmp*diag(lambda_y)*phitmp'+sigma*eye(length(t_y{i})));
                pctmp = repmat(pc_x(i,:),Ky,1);
                pctmp = diag(pctmp(:));
                X{i} = pctmp*repmat(phitmp,1,Kx)';
            end
        elseif method == 2
            for i = 1:n
                phitmp = interp1(out1y',phi_y,t_y{i}','spline');
                W{i} = pinv(phitmp*diag(lambda_y)*phitmp'+sigma*eye(length(t_y{i})));
                phitmp = reshape(repmat(phitmp,Kx,1),length(t_y{i}),Kx*Ky);
                pctmp = repmat(pc_x(i,:),Ky,1)';
                pctmp = diag(pctmp(:));
                X{i} = pctmp*phitmp';
            end
        end
        
    end

        