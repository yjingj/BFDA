%[val] = setVal(X,varname,val)
% This function sets the struct from setOptions() or cell array from FPCA(),
% FPCreg() to the value you want. Those for cell arrays are not recommended,
% but this function would be handy for you to use setOptions() to set 
% most of the default values and use setVal() to change couple of them.
%
%Input: 
%      X:  a struct from setOptions() or cell array from FPCA() or FPCreg()
%varname:  a character string that defines the variable you want from X
%    val:   value to be set to X  with name of "varname"        
%
%Output:
%      X:   updated cell array X
%
%     ex:   px =setOptions();  
%           now, change regular status to be "2" instead of "0"
%           px = setVal(px,'regular',2);
%           
% see also getVal, names
function [X] = setVal(X,varname, val)
      id = strmatch(varname,names(X),'exact');
      if length(id) == 0
          fprintf(1,'Invalid variable name!Use function "names" to get a valid name!\n');
          return;
      else
         if iscell(X) == 1
           X{id} = val;
         elseif isstruct(X) == 1
           X = setfield(X, varname,val);
         else
           fprintf(1, 'Warning: only cell array or struct is supported!\n'); 
           return;
         end
      end
end
