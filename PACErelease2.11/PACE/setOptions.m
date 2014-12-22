%[param] = setOptions(varargin)
% This function sets the optional input arguments for the function PCA().
% =====
% Input: 
% =====
%
% This function allows 3 different ways for input:
% 
% 1)   If no input argument is provided, use default values (regular = 0):
%      ex: 
%      >> p = setOptions() 
%
% 2)   Provide the  name of the argument following by the value to be assign
%      to this argument. The order of the input argument pairs (argument name, value)
%      does not have to follow the same order of the input order like those defined
%      in PCA() or showOptionsNames(). Any undefined input arguments will be followed
%      the same as the default ones.
%
%      ex: 
%      >> showOptionNames()   %see what argument names you can define
%                            %and their orders 
%      >> p = setOptions('regular',2, 'bwmu',3);
%        means set regular = 2, bwmu = 3 and everything else
%        will be using the default values
% 
% 3)   Provide the argument value directly, but ORDER matters. 
%      Any undefined argument values will be set to the default ones.
%   
%      ex:
%      >> showOptionNames()   %see what argument names you can define
%                            %and their orders 
%      >> p = setOptions(3,[],[3 3])
%      means set bwmu = 3, bwmu_gcv =[] and bwxcov = [3 3]
%      everything else will be using the default values
%       
% ======
% Output: a struct array that contains all user-defined options for PCA().
% ======
% 
% See also PCA, showOptionNames, setVal, getVal
function [param] = setOptions(varargin)

        %set default options
        param = struct('bwmu',0,'bwmu_gcv', 1, 'bwxcov',[0 0],...
	               'bwxcov_gcv',1,'ntest1',30,'ngrid1',30,...
	               'selection_k','BIC1', 'FVE_threshold',0.85,...
                   'maxk',20, 'control','auto','regular',[],...
	               'error',1,'ngrid',51,'method','CE',...
      	           'shrink',0,'newdata',[],'kernel','gauss', 'numBins',[],...
                   'yname', [], 'screePlot', 0, 'designPlot',0, 'corrPlot',0,...
                   'rho','cv', 'family','binomial','param',[],'verbose','on','xcov',[]);
        if nargin > 0
            op_names = {'bwmu','bwmu_gcv','bwxcov', ...
	                'bwxcov_gcv','ntest1','ngrid1',...
	                'selection_k','FVE_threshold', 'maxk',...
                        'control','regular',...
	                'error','ngrid','method',...
                        'shrink','newdata','kernel', 'numBins', 'yname'...
                        ,'screePlot','designPlot', 'corrPlot', 'rho','family','param','verbose','xcov'};
            paramLen = length(struct2cell(param)); 
            if ischar(varargin{1})  %match and set arguments by names
                if nargin > 2*paramLen
                      fprintf(1,'Warning: Too many input arguments!Reset to default values!\n');     
                elseif mod(nargin, 2) == 0
                      for i = 1:2:nargin
	                     id = strmatch(varargin{i},op_names,'exact');
                         if isempty(id)
                           fprintf(1,['Warning: "' varargin{i} '" is an invalid name and it will be ignored!Type showOptionNames() for more details.\n']);
                         else
		                   param = setVal(param, op_names{id}, varargin{i+1});   
                         end
                      end
                else
		           fprintf(1,'Warning:Either option name or its corresponding value is missing!Reset to default values!\n');
                end
	        else  %set arguments directly
               if nargin <= paramLen
	               
                    for i = 1:nargin
	                  param = setVal(param, op_names{i}, varargin{i});
                    end
               else
	               fprintf(1,'Warning: Too many input arguments!Reset to default values!\n');
               end      
           end

        end

end
