%[op_names] = showDerOptionNames()
%This function display the optional input argument names for PCA()
%It will be helpful if you want to use setOptions() to set the values
%corresponding to these names.
%=======
%Output:
%======
%op_names:    a 19 * 1 cell array contains the names of the input arguments
%             for PCA)()
%
% ex:
% >> showDerOptionNames()
%see also PCder, FPCder, setDerOptions
function [op_names] = showDerOptionNames()
     fprintf(1, '\nHere is a list of optional input argument names for the PCA():\n');

    op_names = {'bwmu','bwmu_gcv','bwxcov', ...
                'bwxcov_gcv','ntest1','ngrid1',...
                'selection_k','FVE_threshold', 'maxk',...
                'control','regular',...
                'error','ngrid','method',...
                'shrink','newdata','kernel', 'npoly', 'nder', 'method_int', 'numBins', 'yname','screePlot','designPlot','rho','verbose'};

end
