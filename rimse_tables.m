load('StatSimMay14.mat');
load('Z_SP_Kernel.mat');

%%
E_signal = NaN(5, 100);
E_mean = NaN(5, 100);
E_cov = NaN(4, 100);

E_signal_NS = NaN(5, 100);
E_mean_NS = NaN(5, 100);
E_cov_NS = NaN(4, 100);


for(i = 1 : 100)
    
    X = XSim(:, :, i);
    Y = YSim(:, :, i);  
    X_NS = XSim_NS(:, :, i);
    Y_NS = YSim_NS(:, :, i);
    
    Zp  = ZpSim(:, :, i) + repmat(mupSim(:, i), 1, 100);
    Zp_NS  = ZpSim_NS(:, :, i) + repmat(mupSim_NS(:, i), 1, 100);

    E_signal(:, i) = [rimse(t, X, ZSim_SP(:, :, i), 0), rimse(t, X, ZSim_Kern(:, :, i), 0),...
                rimse(t, X, ZSim(:, :, i), 0), rimse(t, X, Zp, 0), rimse(t, X, ZSim_OP(:, :, i), 0)];
    E_signal_NS(:, i) = [rimse(t, X_NS, ZSim_SP_NS(:, :, i), 0), rimse(t, X_NS, ZSim_Kern_NS(:, :, i), 0),...
                rimse(t, X_NS, ZSim_NS(:, :, i), 0), rimse(t, X_NS, Zp_NS, 0), rimse(t, X_NS, ZSim_OP_NS(:, :, i), 0)];
     
    E_mean(:, i) = [rimse(t, mut', mean(ZSim_SP(:, :, i), 2), 0), rimse(t, mut', mean(ZSim_Kern(:, :, i), 2), 0), ...
        rimse(t, mut', muSim(:, i), 0),  rimse(t, mut', mupSim(:, i), 0), rimse(t, mut', mean(ZSim_OP(:, :, i), 2), 0)];
    E_mean_NS(:, i) = [rimse(t, mut_NS', mean(ZSim_SP_NS(:, :, i), 2), 0), ...
        rimse(t, mut_NS', mean(ZSim_Kern_NS(:, :, i), 2), 0), ...
        rimse(t, mut_NS', muSim_NS(:, i), 0),  rimse(t, mut_NS', mupSim_NS(:, i), 0),...
        rimse(t, mut_NS', mean(ZSim_OP_NS(:, :, i), 2), 0)];    
    
    E_cov(:, i) = [rimse(t, C, cov(Y'), 2), rimse(t, C, iKSim(:, :, i), 2), ...
                 rimse(t, C, iKpSim(:, :, i), 2), rimse(t, C, cov(ZSim(:, :, i)'), 2)];
    E_cov_NS(:, i) = [rimse(t, C_NS, cov(Y_NS'), 2), rimse(t, C_NS, iKSim_NS(:, :, i), 2), ...
                         rimse(t, C_NS, iKpSim_NS(:, :, i), 2), rimse(t, C_NS, cov(ZSim_NS(:, :, i)'), 2)];  
                            
end

save('rimse_error.mat', 'E_signal', 'E_signal_NS', 'E_mean', 'E_mean_NS', 'E_cov', 'E_cov_NS');

%%
mean(E_signal, 2)
std(E_signal, 0,2) ./ 10

mean(E_signal_NS, 2)
std(E_signal_NS, 0, 2)./ 10

mean(E_mean, 2)
std(E_mean, 0, 2)./ 10

mean(E_mean_NS, 2)
std(E_mean_NS, 0, 2)./ 10

mean(E_cov, 2)
std(E_cov, 0, 2)./ 10

mean(E_cov_NS, 2)
std(E_cov_NS, 0, 2)./ 10

%%







