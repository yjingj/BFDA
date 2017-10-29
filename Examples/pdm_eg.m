
%% Example of testing the goodness-of-fit 
% by the PDM method proposed by Yuan and Johnson (2012)

addpath(genpath(cat(2, pwd, '/BFDA_SourceCode'))) % replace pwd by the directory of your MATLAB package

%% Simulate a matrix of residuals
w = normrnd(0, 1, 1000, 1, 1000);

% hist(reshape(w(1, 1, :), 1, 1000))

% Call the test of goodness-of-fit
[pmin, pdm, df] = pdm_test(w, 4, 100);

display('test p-value:');
display(pmin);
