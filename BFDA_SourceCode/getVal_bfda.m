%[val] = getVal_bfda(X,varname)
%This function returns the variable from an object returned by BFDA
%based on the name of the variable
%Input: 
%      X:  a cell array that is returned by BFDA
%varname:  a character string that defines the variable you want from X
%Output:
%   val: variable from X with name of "varname"        
% ex:   no_opt = getVal_bfda(X,'rho'); returns estimate for rho
%       or x =setOptions_bfda();  bwmu = getVal_bfda(x, 'bwmu');

function [val] = getVal_bfda(X,varname)
   
      id = strncmp(varname,names_bfda(X),9);
      if sum(id)==0
          fprintf(1,'Invalid variable name!Use function "names" to get a valid name!\n');
          val = [];
      else
         if iscell(X) == 1 && sum(id) == 1
           val = X{id};
         elseif isstruct(X) == 1 && sum(id) == 1
           val = X.(varname);
         elseif sum(id) > 1
           fprintf(1, 'Warning: multiple matches for field name!\n');
         else
           fprintf(1, 'Warning: only cell array or struct is supported!\n');
	   val = [];
         end
      end
end
