% This function sets the optional input arguments for the functions bhmcmc(), babf_mcmc().
% =====
% Input: 
% =====
%
% This function allows 3 different ways for input:
% 
% 1)   If no input argument is provided, use default values (mat = 1):
%      ex: 
%      >> param = setOptions_bfda() 
%
% 2)   Provide the  name of the argument following by the value to be assign
%      to this argument. The order of the input argument pairs (argument name, value)
%      does not have to follow the same order of the input order like those defined
%      in BFDA() or showOptionsNames_bfda(). Any undefined input arguments will be followed
%      the same as the default ones.
%
%      ex: 
%      >> showOptionNames_bfda()   %see what argument names you can define
%                            %and their orders 
%      >> param = setOptions_bfda('regular',2, 'bwmu',3);
%        means set regular = 2, bwmu = 3 and everything else
%        will be using the default values
% 
% 3)   Provide the argument value directly, but ORDER matters. 
%      Any undefined argument values will be set to the default ones.
%   
%      ex:
%      >> showOptionNames_bfda()   %see what argument names you can define
%                            %and their orders 
%      >> param = setOptions_bfda('bhm',[],1)
%      means set smethod = 'bhm', cgrid = [], and mat = 1
%      everything else will be using the default values
%       
% ======
% Output: a struct array that contains all user-defined options 
%           for functions bhmcmc(), babf_mcmc() ...
% ======
% 
% See also showOptionNames_bfda, setVal_bfda, getVal_bfda

function [param] = setOptions_bfda(varargin)

        %set default options
        param = struct('smethod', 'babf', 'Burnin', 2000, 'M', 10000, 'cgrid', 1, 'mat', 1, ...
            'Sigma_est', [], 'mu_est', [], 'nu', [],  ...
            'delta', 5, 'c', 1, 'w', 1, 'ws', 0.1, 'pace', 1, ...
            'tau', [], 'Regress', 0 , 'm', 20, 'eval_grid', [], ...
            'a', 0.001, 'b', 0.001, 'lamb_min', 0.9, ...
            'lamb_max', 0.99, 'lamb_step', 0.01, 'resid_thin', 10, 'tol', 0);
        
        if nargin > 0
            op_names = {'smethod', 'Burnin', 'M', 'cgrid', 'mat',  ...
            'Sigma_est', 'mu_est', 'nu',  ...
            'delta',  'c', 'w', 'ws', 'pace', ...
            'tau', 'Regress', 'm', 'eval_grid', ...
            'a', 'b', 'lamb_min', 'lamb_max', 'lamb_step', 'resid_thin', 'tol'};
        
            paramLen = length(struct2cell(param)); 
            if ischar(varargin{1})  %match and set arguments by names
                if nargin > 2*paramLen
                      fprintf(1,'Warning: Too many input arguments!Reset to default values!\n');     
                elseif mod(nargin, 2) == 0
                      for i = 1:2:nargin
	                     id = strncmp(varargin{i}, op_names, 9);
                         if sum(id) == 0
                           fprintf(1,['Warning: "' varargin{i} '" is an invalid name and it will be ignored!Type showOptionNames_bfda() for more details.\n']);
                         elseif sum(id) == 1
		                   param = setVal_bfda(param, op_names{id}, varargin{i+1}); 
                         else
                            fprintf(1,['Warning: Find multiple match names for "' varargin{i} '"! Type showOptionNames_bfda() for more details.\n']); 
                         end
                      end
                else
		           fprintf(1,'Warning:Either option name or its corresponding value is missing!Reset to default values!\n');
                end
	        else  %set arguments directly
               if nargin <= paramLen
	               
                    for i = 1:nargin
	                  param = setVal_bfda(param, op_names{i}, varargin{i});
                    end
               else
	               fprintf(1,'Warning: Too many input arguments!Reset to default values!\n');
               end      
           end

        end

end
