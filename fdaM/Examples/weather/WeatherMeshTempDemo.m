addpath('/Users/jamesramsay/Documents/MATLAB/fdaM/examples/weather')
addpath('/Users/jamesramsay/Documents/MATLAB/fdaM')

load WeatherMesh

% fhdl = figure(1);
% 
% plot(basisobj)
% set(fhdl, 'Position', [000 518 800 800])
fhdl = figure(2);
set(fhdl, 'Position', [200 518 1000 800])
for i=weekindex
    plot_FEM(weektempfd(i), AZ, EL)
    axis([-150,-50,40,90,-45,25])
    xlabel('\fontsize{24} Longitude')
    ylabel('\fontsize{24} Latitude')
    view(AZ,EL)
    caxis([-45,25])
    shading interp
    title(['\fontsize{24} Temperature for Week ',num2str(i),' '])
    pause(1);
end
