function slidemat = slides(x, breaks, rates, nderiv)
%  SLIDE evaluates a slide basis
%  Slide defines a series of exponential decays with rates defined by
%  values in argument X over intervals defined by BREAKS

%  Last modified 14 March 2010

x = x(:);

sizebreaks = size(breaks);
if sizebreaks(1) > 1 && sizebreaks(2) > 1
    error('Argument BREAKS is not a vector.');
end
breaks = breaks(:);

%  check break point values

nbreaks = length(breaks);
if nbreaks < 2, error('Number of breaks less than 2.'); end

%  check BREAKS for being increasing

if min(diff(breaks)) < 0 
    breaks = sort(breaks); 
end

if min(diff(breaks)) == 0
    error('Values are not unique.');
end

%  check rates

nrates = length(rates);
if nrates ~= nbreaks - 1
    error('Number decay rates not equal to NBREAKS - 1.');
end

%  set up an order 1 basis

rangeval = [breaks(1), breaks(nbreaks)];

StepBasis = create_bspline_basis(rangeval, nbreaks-1, 1, breaks);

slidemat = eval_basis(x, StepBasis);
for ibreak=2:nbreaks
    index = find(breaks(ibreak-1) <= x & x < breaks(ibreak));
    if ~isempty(index)
        ratei = rates(ibreak-1);
        slidei = ...
            (-ratei).^nderiv.*exp(-ratei.*(x(index)-breaks(ibreak-1)));
        slidemat(index,ibreak-1) = slidemat(index,ibreak-1).*slidei;
    end
end


