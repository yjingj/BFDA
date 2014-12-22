%[val] = getVal(X,varname)
%This function returns the variable from an object returned by FPCA
%based on the name of the variable
%Input: 
%      X:  a cell array that is returned by FPCA
%varname:  a character string that defines the variable you want from X
%Output:
%   val: variable from X with name of "varname"        
% ex:   no_opt = getVal(X,'no_opt'); returns the number of optimal principal components
%       or x =setOptions();  bwmu = getVal(x, 'bwmu');
function [val] = getVal(X,varname)
   
      id = strmatch(varname,names(X),'exact');
      if isempty(id)
          fprintf(1,'Invalid variable name!Use function "names" to get a valid name!\n');
          val = [];
      else
         if iscell(X) == 1
           val = X{id};
         elseif isstruct(X) == 1
           val = getfield(X, varname);
         else
           fprintf(1, 'Warning: only cell array or struct is supported!\n');
	   val = [];
         end
      end
end
