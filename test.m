%%
GausData = gaussian_sim(t, n, sf, snr, nu, rho, dense);

t = GausData.t;
C = GausData.C;
D = GausData.D;
X0 = GausData.X0;
Y0 = GausData.Y0;
mut = GausData.mut;

X = GausData.X;
Y = GausData.Y;
T = GausData.T;

X_cgrid = GausData.X_cgrid;
Y_cgrid = GausData.Y_cgrid;
T_cgrid = GausData.T_cgrid;

BHM_cgrid = bhmcmc(Y_cgrid, T_cgrid, delta, 2, ...
        5000, 5000, mat, Sigma_est, mu_cgrid, [], 0.5, 1); 
    rmse(BHM_cgrid.iK, C) 
    rmse(Sigma_est, C)
    rmse(COR(BHM_cgrid.iK), COR(C))
    rmse(COR(Sigma_est), COR(C))
    rmse(BHM_cgrid.Z, X0)
   