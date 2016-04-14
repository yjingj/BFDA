fdaMpath = 'D:\spencerg\statmtds\splines\fda\RForge\fda\pkg\inst\Matlab\fdaM' ; 
addpath(fdaMpath) ; 

addpath('c:/PROGRA~1/r/R-27~1.0/library/R.matlab/externals');


fdaMdirs = {'examples', 'examples/gait', 'examples/gdp', ... 
    'examples/goodsindex', 'examples/growth', 'examples/handwrit', ... 
    'examples/lip', 'examples/mathtest', 'examples/melanoma', ...
    'examples/pinch', 'examples/refinery', 'examples/weather', ...
    'ftpsite', 'ftpsite/temp'};

for id = 1: length(fdaMdirs) 
    ifp = fullfile(fdaMpath, fdaMdirs{id}) ; 
    addpath(ifp)  ;  
end 
