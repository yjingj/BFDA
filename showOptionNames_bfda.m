%[op_names] = showOptionNames_bfda()
%This function display the optional input argument names for BFDA()
%It will be helpful if you want to use setOptions_bfda() to set the values
%corresponding to these names.
%=======
%Output:
%======
%op_names:    a 1 * 14 cell array contains the names of the input arguments
%             for BFDA()
%
% ex:
% >> showOptionNames_bfda()
%see also setOptions_bfda

function [op_names] = showOptionNames_bfda()
     fprintf(1, '\nHere is a list of optional input argument names for the BFDA():\n');
     op_names = {'smethod', 'Burnin', 'M', 'cgrid', 'mat',  ...
            'Sigma_est', 'mu_est', 'nu',  ...
            'delta',  'c', 'w', 'ws', 'pace', 'tau', ...
            'Regress', 'm', 'eval_grid', 'a', 'b', ...
            'lamb_min', 'lamb_max', 'lamb_step'};
end
