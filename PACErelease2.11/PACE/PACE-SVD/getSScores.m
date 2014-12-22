function [xi_est, yi_est, ss_var, rho_opt, sig1]=getSScores(res, x, t_x, y, t_y, nsvd, error, method, shrink,  regular, rho)

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
    rho_opt=[];        

    n = length(x);
    sig1 = [sigma_x sigma_y];

    %Compute iterative residuals
    if rho ~= -1
      for j = 1:2
        SE_i = inf*ones(n,2);
        %get fitted curves based on t{i} for a given TOL(j)
        xo = getOriCurves(x, t_x, mu_x, sc_x, lambda_x, sigma_x, sig1(1), nsvd, error(1), method, shrink, out1x, regular);        
        yo = getOriCurves(y, t_y, mu_y, sc_y, lambda_y, sigma_y, sig1(2), nsvd, error(2), method, shrink, out1y, regular);        
        for i = 1:n
           SE_i(i,1) = mean((x{i}-xo{i}).^2);
           SE_i(i,2) = mean((y{i}-yo{i}).^2);
        end
        sig1 = mean(SE_i);
      end
    end

    if isnumeric(rho) 
       if rho >= 0 || rho == -1
         rho_opt = rho;
       else
	 fprintf(1,'Warning: rho should not be negative! Reset it to cv choice now!\n');
         rho = 'cv';
       end
    end

    if ischar(rho) && strcmp('cv',rho)

        %Compute gamma
        %fprintf(1, ['The new sigma after 2 iterations is ' num2str(sig1) '\n']);
        T_x = range(cell2mat(t_x));
        gamma_x = ((trapz(out1x, mu_x.^2)+sum(lambda_x))/T_x)^(0.5);
        T_y = range(cell2mat(t_y));
        gamma_y = ((trapz(out1y, mu_y.^2)+sum(lambda_y))/T_y)^(0.5);
        %fprintf(1,['Gamma = ' num2str(gamma) '\n']);
        %alpha = linspace(0.001*gamma, 0.05*gamma, 50);
        %alpha = linspace(0.005,0.22,50);
        alpha = linspace(0.01,0.44,50);
        rho = [];
        rho(:,1) = gamma_x*alpha;
        rho(:,2) = gamma_y*alpha;

        ni = zeros(n,2);
        tjID = zeros(n,2);
        for i = 1:n
            for k=1:2
                ni(i,1) = length(t_x{i});
                ni(i,2) = length(t_y{i});
                if ni(i,k) > 1
                    tjID(i,k) = mysample(1:ni(i,k),1,0); %find random index of random time point to be used in the CV prediction part
                    %for ni >= 2
                end
            end
        end

        %find optimal rho from subjects with ni >= 2
        if all(ni(:,1) == 1) || any(error == 0) || all(ni(:,2) == 1) 
           rho_opt = min(rho);
        else
           rho_opt = cv_srho(res, x, t_x, y, t_y, nsvd, sig1, method, shrink, regular, rho, ni, tjID);
        end
    end

    [xi_est, yi_est, ss_var]=getSScores1(res,  x, t_x, y, t_y, nsvd, sig1, error, method, shrink, regular, rho_opt);


end
