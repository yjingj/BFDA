function [xi_est,  yi_est, ss_var]=getSScores1(res,  x, t_x, y, t_y, nosvd, sigma_new, error, method, shrink, regular, rho)

mu_x=getVal(getVal(res,'xx'),'mucopy');
sc_x=getVal(res,'sc_x');
lambda_x = getVal(res,'lambda_x');
sigma_x = getVal(getVal(res,'xx'),'sigma');
out1x = getVal(res,'out_x');
mu_y=getVal(getVal(res,'yy'),'mucopy');
sc_y=getVal(res,'sc_y');
lambda_y = getVal(res,'lambda_y');
sigma_y = getVal(getVal(res,'yy'),'sigma');
out1y = getVal(res,'out_y');
lambda = getVal(res, 'lambda');
[muxSub, scxSub] = convertMuPhi(t_x, out1x, mu_x, sc_x, regular);
[muySub, scySub] = convertMuPhi(t_y, out1y, mu_y, sc_y, regular);
ncohort = length(x);

if strcmp(method,'IN')
    ss_var = [];
else
    ss_var = cell(1,ncohort);
end
%fprintf(1,'Start calculating the PC scores \n');
%update \xi

if error==1

    sigma1 = max(sigma_new, rho);

    if regular == 2 && strcmp(method,'CE')

        xx= reshape(cell2mat(x), length(x{1}), ncohort)';
        yy= reshape(cell2mat(y), length(y{1}), ncohort)';

        A22 = scxSub*diag(lambda_x)*scxSub'+sigma1(1)*eye(length(x{1}));
        B22 = scxSub*diag(lambda)*scySub';
        D22 = scySub*diag(lambda_y)*scySub'+sigma1(2)*eye(length(y{1}));
        A22i = pinv(A22-B22*pinv(D22)*B22');
        D22i = pinv(D22-B22'*pinv(A22)*B22);
        B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
        sig22i = [A22i B22i; B22i' D22i];
        sig12 = [diag(lambda_x)*scxSub' diag(lambda)*scySub'; diag(lambda)*scxSub' diag(lambda_y)*scySub'];
        sig11 = [diag(lambda_x) diag(lambda); diag(lambda) diag(lambda_y)];
        mu_ss = (sig12 * sig22i * ([xx yy] - repmat([muxSub muySub],ncohort,1))')';
        xi_est = mu_ss(:,1:nosvd);
        yi_est = mu_ss(:,(nosvd+1):end);

        %         y_predOrig = MU+xi_est*phiSub';
        %         y_predOrig = num2cell(y_predOrig,2);
        C = sig11 - sig12 * sig22i * sig12';
        for i = 1:ncohort
            ss_var{i} = C;
        end

    else
%         x_predOrig = cell(1,ncohort);
%         y_predOrig = cell(1,ncohort);
        xi_est = zeros(ncohort, nosvd);
        yi_est = zeros(ncohort, nosvd);
        zeta_xest = xi_est;
        zeta_yest = yi_est;
        scxi= scxSub;
        mux_i = muxSub;
        scyi= scySub;
        muy_i = muySub;
        for i = 1:ncohort

            if regular ~= 2
                scxi = scxSub{i};
                mux_i = muxSub{i};
                scyi = scySub{i};
                muy_i = muySub{i};
            end
            xi= x{i};
            yi= y{i};
            if strcmp(method,'CE')

                A22 = scxi*diag(lambda_x)*scxi'+sigma1(1)*eye(length(xi));
                B22 = scxi*diag(lambda)*scyi';
                D22 = scyi*diag(lambda_y)*scyi'+sigma1(2)*eye(length(yi));
                A22i = pinv(A22-B22*pinv(D22)*B22');
                D22i = pinv(D22-B22'*pinv(A22)*B22);
                B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
                sig22i = [A22i B22i; B22i' D22i];
                sig12 = [diag(lambda_x)*scxi' diag(lambda)*scyi'; diag(lambda)*scxi' diag(lambda_y)*scyi'];
                sig11 = [diag(lambda_x) diag(lambda); diag(lambda) diag(lambda_y)];
                mu_ss = sig12 * sig22i * ([xi yi] - [mux_i muy_i])';
                xi_est(i,:) = mu_ss(1:nosvd);
                yi_est(i,:) = mu_ss((nosvd+1):end);

                ss_var{i} = sig11 - sig12 * sig22i * sig12';
            elseif strcmp(method,'IN')
                mx=length(xi);
                my=length(yi);
                for k=1:nosvd
                    prodx=(xi-mux_i).*scxi(:,k)';
                    prody=(yi-muy_i).*scyi(:,k)';
                    if shrink == 0
                        %xi_est(i,k) = romb(t{i},prod);
                        xi_est(i,k) = trapz(t_x{i},prodx);
                        yi_est(i,k) = trapz(t_y{i},prody);
                        %zeta_est(i,k)=trapzoid(prod,t{i});
                    else
                        %zeta_est(i,k) = romb(t{i},prod);
                        zeta_xest(i,k) = trapz(t_x{i},prodx);
                        zeta_yest(i,k) = trapz(t_y{i},prody);
                        xi_est(i,k)=lambda_x(k)*zeta_xest(i,k)/(lambda_x(k)+sigma_x/mx);
                        yi_est(i,k)=lambda_y(k)*zeta_yest(i,k)/(lambda_y(k)+sigma_y/my);
                    end
                end
            end
%             y_predOrig{i} = mu_i+xi_est(i,:)*phii';
        end
    end
elseif error==0
    if regular == 2 && strcmp(method ,'CE')
        xx= reshape(cell2mat(x), length(x{1}), ncohort)';
        yy= reshape(cell2mat(y), length(y{1}), ncohort)';

        A22 = scxSub*diag(lambda_x)*scxSub';
        B22 = scxSub*diag(lambda)*scySub';
        D22 = scySub*diag(lambda_y)*scySub';
        A22i = pinv(A22-B22*pinv(D22)*B22');
        D22i = pinv(D22-B22'*pinv(A22)*B22);
        B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
        sig22i = [A22i B22i; B22i' D22i];
        sig12 = [diag(lambda_x)*scxSub' diag(lambda)*scySub'; diag(lambda)*scxSub' diag(lambda_y)*scySub'];
        sig11 = [diag(lambda_x) diag(lambda); diag(lambda) diag(lambda_y)];
        mu_ss = (sig12 * sig22i * ([xx yy] - repmat([muxSub muySub],ncohort,1))')';
        xi_est = mu_ss(:,1:nosvd);
        yi_est = mu_ss(:,(nosvd+1):end);

%         y_predOrig = MU+xi_est*phiSub';
%         y_predOrig = num2cell(y_predOrig,2);
        C = sig11 - sig12 * sig22i * sig12';
        for i = 1:ncohort
            ss_var{i} = C;
        end
    else
%         y_predOrig = cell(1,ncohort);
        xi_est = zeros(ncohort, nosvd);
        yi_est = zeros(ncohort, nosvd);
        scxi= scxSub;
        scyi= scySub;
        mux_i = muxSub;
        muy_i = muySub;
        for i = 1:ncohort
            if regular ~= 2
                scxi = scxSub{i};
                scyi = scySub{i};
                mux_i = muxSub{i};
                muy_i = muySub{i};
            end
            xi= x{i};
            yi= y{i};
            if strcmp(method,'CE')
                A22 = scxi*diag(lambda_x)*scxi';
                B22 = scxi*diag(lambda)*scyi';
                D22 = scyi*diag(lambda_y)*scyi';
                A22i = pinv(A22-B22*pinv(D22)*B22');
                D22i = pinv(D22-B22'*pinv(A22)*B22);
                B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
                sig22i = [A22i B22i; B22i' D22i];
                sig12 = [diag(lambda_x)*scxi' diag(lambda)*scyi'; diag(lambda)*scxi' diag(lambda_y)*scyi'];
                sig11 = [diag(lambda_x) diag(lambda); diag(lambda) diag(lambda_y)];
                mu_ss = sig12 * sig22i * ([xi yi] - [mux_i muy_i])';
                xi_est(i,:) = mu_ss(1:nosvd);
                yi_est(i,:) = mu_ss((nosvd+1):end);

                ss_var{i} = sig11 - sig12 * sig22i * sig12';
                
            elseif strcmp(method,'IN')
                for k=1:nosvd
                    prodx=(xi-mux_i).*scxi(:,k)';
                    prody=(yi-muy_i).*scyi(:,k)';
                    %xi_est(i,k)=trapzoid(prod,t{i});
                    %xi_est(i,k) = romb(t{i},prod);
                    xi_est(i,k) = trapz(t_x{i},prodx);
                    yi_est(i,k) = trapz(t_y{i},prody);
                end
            end
%             y_predOrig{i} = mu_i+xi_est(i,:)*phii';
        end
    end
end
end
