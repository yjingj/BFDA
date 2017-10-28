%[val] = setVal_bfda(X,varname,val)
% This function sets the struct from setOptions_bfda() or cell array from BFDA(). 
% Those for cell arrays are not recommended,
% but this function would be handy for you to use setOptions() to set 
% most of the default values and use setVal() to change couple of them.
%
%Input: 
%      X:  a struct from setOptions_bfda() or cell array from BFDA()
%varname:  a character string that defines the variable you want from X
%    val:   value to be set to X  with name of "varname"        
%
%Output:
%      X:   updated cell array X
%
%     ex:   px =setOptions_bfda();  
%           now, change babf status to be "1" instead of "0"
%           px = setVal_bfda(px,'babf',1);
%           
% see also getVal_bfda, names_bfda

function [X] = setVal_bfda(X,varname, val)
      id = strncmp(varname,names_bfda(X),9);
      if sum(id) == 0
          fprintf(1,'Invalid variable name!Use function "names" to get a valid name!\n');
          return;
      else
         if iscell(X) == 1 && sum(id) == 1
           X{id} = val;
         elseif isstruct(X) == 1 && sum(id) == 1
             X.(varname) = val;
           %X = setfield(X, varname,val);
         else
           fprintf(1, 'Warning: only cell array or struct is supported!\n'); 
           return;
         end
      end
end
