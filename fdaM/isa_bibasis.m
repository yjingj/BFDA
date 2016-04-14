function isabasis = isa_bibasis(bibasisobj)
%  ISA_BIBASIS  checks a struct object for fields for bibasis objects

%  last modified 13 December 2012

isabasis = 1;
if ~strcmp(class(bibasisobj),'bibasis')
    isabasis = 0;
end

