%[n] = names_bfda(X)
%This function returns the names of the object returned by BFDA.m
%Input: X is a cell array with a name cell array associating with it
%       or a struct, e.g. object from BFDA(), object from setOptions_bfda() etc.
%Output: n is a cell array of names for X
%ex: names_bfda(X) 

function [n] = names_bfda(X)
     if iscell(X) == 1
        n = X{end};
     elseif isstruct(X) == 1
        n = fieldnames(X)';
     else
        fprintf(1,'Warning: only cell array or struct is supported!\n');
        n = [];
     end
end
