function plot_FEM(FEMfd, AZ, EL, clims, label)
% PLOT  Plots a FEM object FEMFD over its triangular mesh.
% The coefficient vector for FEMfd object defines the heights of the 
% surface at the vertices.  The vertices are in
% NBASIS by 2 matrix P extracted from the FEM basis object.
%
% Plotting is done using Matlab function trisurf.
%
% Last modified on 2 November 2014

if nargin < 3, EL =  30.0;  end
if nargin < 2, AZ = -37.5;  end

coefmat = full(getcoef(FEMfd));
nsurf   = size(coefmat,2);

basisobj = getbasis(FEMfd);

params = getbasispar(basisobj);

p = params.p;
t = params.t(:,1:3);

for isurf=1:nsurf
    trisurf(t, p(:,1), p(:,2), coefmat(:,isurf))
    view(AZ,EL)
    colorbar
    if nargin >=4
        caxis(clims)
    end
    shading interp
    xlabel('\fontsize{13} X')
    ylabel('\fontsize{13} Y');
    if nsurf > 1
        if nargin < 5
            title(['\fontsize{16} Sample ',num2str(isurf)]);
        else
            title(['\fontsize{16} ',label(isurf,:)]);
        end
        pause
    end
end


