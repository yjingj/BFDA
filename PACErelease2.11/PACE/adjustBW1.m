function [bopt] = adjustBW1(kernel,bopt, npoly, nder, regular,verbose)

  if strcmp(kernel,'gauss')
        if regular == 2
           bwmu_fac = [1.1 0.8 0.8];
        else
           bwmu_fac = [1.1 1.2 2];
        end
        if nder > 2
	  facID = 3;
        elseif nder >= 0 && nder <= 2
          facID = nder+1;
        else
	  facID = 1;
        end
        bopt = bopt*bwmu_fac(facID);
        if strcmp(verbose, 'on') == 1
          fprintf(1,['Adjusted GCV bandwidth choice for mean function (npoly = ' num2str(npoly) '): ' num2str(bopt) '\n']);
        end
  end

  if strcmp(kernel,'epan')
        if regular == 2
           bwmu_fac = [1.1 1.0 1.1];
        else
           bwmu_fac = [1.1 1.2 1.5];
        end
        if nder > 2
	  facID = 3;
        elseif nder >= 0 && nder <= 2
          facID = nder+1;
        else
	  facID = 1;
        end
        bopt = bopt*bwmu_fac(facID);
        if strcmp(verbose, 'on') == 1
          fprintf(1,['Adjusted GCV bandwidth choice for mean function (npoly = ' num2str(npoly) '): ' num2str(bopt) '\n']);
        end
  end

end
