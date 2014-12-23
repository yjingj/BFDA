SmoothingMatlabCode
===================

## Matlab code to smoothing functional data

### Grid approach related functions
- **bhmcmc.m**
Main function to implement MCMC

- **SimulationScript.m**
Example simulation studies

- **CaseStudyScript.m**
Multiple case study script

- **plot...m**
Script to generate plots

### Basis function approach related code

### Regression related code

### Simulation related functions
- **gaussian_sim.m**
Function to generate simulated stationary functional data with Matern covariance.

- **gaussian_simNS.m**
Function to generate non-stationary functional data with covariance trasformed from a Matern class.

### Other function used in the code
- **Matern.m**
Function to calculate Matern covariance given distance matrix.

- **COR.m**
Function to output correlation matrix given covariance matrix.

- **rmse.m**
Function to calculate RMSE.

- **rimse.m**
Function to calculate RIMSE.

- **mychol.m**
Function to output a matrix L such that LL' = input covariance matrix.

- ** myiwishrnd.m **
Function to generate sample matrix from a Inverse-Wishart distribution with given degrees of freedom, scale matrix, and dimension.

### Other Matlab Libraries used in the code

- **PACErelease2.11/**
PACE matlab package written by Yao et. al.

- **mcmcdiag/**
Matlab functions to diagonize MCMC convergence.
