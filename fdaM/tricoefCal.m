function tricoef = tricoefCal(p, t)
%  TRICOEFCAL compute the coefficient matrix TRICOEF
%  required to test of a point is indside a triangle

ntri   = size(t,1);

%  compute coefficients for computing barycentric coordinates if
%  needed

tricoef = zeros(ntri,4);
tricoef(:,1) = p(t(:,1),1)-p(t(:,3),1);
tricoef(:,2) = p(t(:,2),1)-p(t(:,3),1);
tricoef(:,3) = p(t(:,1),2)-p(t(:,3),2);
tricoef(:,4) = p(t(:,2),2)-p(t(:,3),2);
detT = tricoef(:,1).*tricoef(:,4) - tricoef(:,2).*tricoef(:,3);
tricoef = tricoef./(detT*ones(1,4));
