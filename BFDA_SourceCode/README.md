
BFDA Matlab Toolbox
===================

## Matlab toolbox for Bayesian functional data analysis by a hierarchical model

### Required Existing MATLAB Libraries for using this toolbox (available under the directory BFDA_SRC/)

- **bspline/**
Matlab toolbox for implementing B-splines, can be downloaded from http://www.mathworks.com/matlabcentral/fileexchange/27374-b-splines.

- **fdaM/**
Matlab toolbox for standard functional data analysis, can be downloaded from http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/.

- **mcmcdiag/**
Matlab toolbox for diagonizing MCMC convergence, can be downloaded from http://becs.aalto.fi/en/research/bayes/mcmcdiag/.

- **PACErelease2.11/**
Matlab toolbox for PACE, written by Yao et. al, can be downloaded from http://www.stat.ucdavis.edu/PACE/download.html.

### Add the path of required MATLAB libraries 
For example, run the following script lines (replace pwd with the directory of your MATLAB libraries)
```
	addpath(genpath(cat(2, pwd, '/bspline')))
	addpath(genpath(cat(2, pwd, '/fdaM'))) 
	addpath(genpath(cat(2, pwd, '/mcmcdiag'))) 
	addpath(genpath(cat(2, pwd, '/PACErelease2.11')))  
```

### Main function of BFDA
- **BFDA_SRC/BFDA.m**
-- Main function for smoothing and mean-covariance estimation.

- **Examples/BayesianSmooth_eg.m**
Example script calling BFDA().

- **Examples/FdRegress_eg.m**
Example script calling fdaM regression functions with output from BFDA().

### Simulation related functions
- **BFDA_SRC/sim_gfd.m**
Function to generate functional data with common/partial grids, stationary/nonstationary covariance.

- **BFDA_SRC/sim_gfd_rgrid.m**
Function to generate functional data with random grids, stationary/nonstationary covariance.



