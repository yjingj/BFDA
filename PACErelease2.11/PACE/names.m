%[n] = names(X)
%This function returns the names of the object returned by FPCA.m
%Input: X is a cell array with a name cell array associating with it
%       or a struct, e.g. object from FPCA(), object from setOptions() etc.
%Output: n is a cell array of names for X
%ex: names(X) 
function [n] = names(X)
     if iscell(X) == 1
        n = X{end};
     elseif isstruct(X) == 1
        n = fieldnames(X)';
     else
        fprintf(1,'Warning: only cell array or struct is supported!\n');
        n = [];
     end
end
