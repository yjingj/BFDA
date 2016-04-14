type = getbasistype(bibasisobj)

%  retrieve the two separate basis objects from a bivariate basis object

%  Last modified 23 November 2012

if ~strcmp(class(bibasisobj), 'bibasis')
    error('BIBASISOBJ is not a bibasis object.');
end

type = bibasisobj.type;
