function [xi_est, xi_var,y_predOrig, rho_opt, sig1]=getScores(y, t, mu, phi, lambda, sigma, noeig, error, method, shrink, out1, regular, rho,verbose)

            
    n = length(y);
    sig1 = sigma;

    %Compute iterative residuals
    if rho ~= -1
      for j = 1:2
        SE_i = inf*ones(1,n);
        %get fitted curves based on t{i} for a given TOL(j)
	yo = getOriCurves(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular);        
        for i = 1:n
           SE_i(i) = mean((y{i}-yo{i}).^2);
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

%    if ischar(rho) && strcmp('cv',rho)
    if ischar(rho) && isempty(strfind(rho, 'cv')) == 0

        %Compute gamma
        %fprintf(1, ['The new sigma after 2 iterations is ' num2str(sig1) '\n']);
        T = range(cell2mat(t));
        gamma = ((trapz(out1, mu.^2)+sum(lambda))/T)^(0.5);
        %fprintf(1,['Gamma = ' num2str(gamma) '\n']);
        %alpha = linspace(0.001*gamma, 0.05*gamma, 50);
        %alpha = linspace(0.005,0.22,50);
        alpha = linspace(0.01,0.22,50);
        rho = gamma*alpha;

        %ni = zeros(1,n);
        %tjID = zeros(1,n);
        %for i = 1:n
        %   ni(i) = length(t{i});
        %   if ni(i) > 1
        %     tjID(i) = mysample(1:ni(i),1,0); %find random index of random time point to be used in the CV prediction part
        %                                      %for ni >= 2
        %   end
        %end
        
        %default is 'cv', non-randomized leave-one-measurement-out
        isRandom = 0;
        if strcmp('cv-random', rho)
	  isRandom = 1;  % randomized leave-one-measurement-out
        end

        [tjID ni] = getTimeID(t,n, isRandom);

        %find optimal rho from subjects with ni >= 2
        if all(ni == 1) || error == 0
           rho_opt = min(rho);
        else
	   rho_opt = cv_rho(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, rho, ni, tjID,verbose);
        end
    end

    [xi_est, xi_var,y_predOrig]=getScores1(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, rho_opt);    


end
  
function [tjID ni] = getTimeID(t, n, isRandom) 

  ni = zeros(1,n);
  tjID = zeros(1,n);

  if isRandom == 1    
     for i = 1:n
        ni(i) = length(t{i});
        if ni(i) > 1
          tjID(i) = mysample(1:ni(i),1,0); %find random index of random time point to be used in the CV prediction part
                                           %for ni >= 2
        end
     end
  else
    
    for i = 1:n
       ni(i) = length(t{i});
       tjID(i) = mod(1000+i,ni(i))+1;      %find non-randomized index of time points
    end
       
  end
  

end

