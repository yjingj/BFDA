BFDA Matlab Toolbox
===================

## Matlab toolbox for smoothing, estiamtion, and regression of functional data by a Bayesian hierarchial model

### Main function to call
- **BFDA.m**

-- Smoothing and mean-covariance estimation

-- Regression with single functional observation

-- Regression with a batch of functional observations

- **SimulationScript.m**
Example simulation studies


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
