BFDA Matlab Toolbox
===================

## Matlab toolbox for Bayesian functional data analysis by a hierarchial model

### Main function to call
- **BFDA.m**
-- Main function for smoothing and mean-covariance estimation.

- **Examples/BayesianSmooth_eg.m**
Example script calling BFDA().

- **Examples/FdRegress_eg.m**
Example script calling fdaM regression functions with output from BFDA().

### Simulation related functions
- **sim_gfd.m**
Function to generate functional data with common/partial grids, stationary/nonstationary covariance.

- **sim_gfd_rgrid.m**
Function to generate functional data with random grids, stationary/nonstationary covariance.

### Other Matlab Libraries used in the code
- **PACErelease2.11/**
Matlab toolbox for PACE, written by Yao et. al.

- **mcmcdiag/**
Matlab toolbox for diagonizing MCMC convergence.

- **bspline/**
Matlab toolbox for implementing B-splines.

- **fdaM/**
Matlab toolbox for standard functional data analysis.

