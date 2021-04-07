clear, clc, close all
addpath(genpath(pwd));
n = 2;
xE = [0.2500    ,0.7500    ,0.5000;
      0.2500   , 0.5000  ,  0.87500];
% xE = [0.2500    ,0.7500    ,0.5000];
bounds = [zeros(n, 1), ones(n, 1)];
surrogate_ = 'constant'; % or use 'constant' to switch to constant surrogate model
func_eval = @(x) -sum((2*x) .* sin(sqrt(500*abs(x))));
xmin = 0.8419 * ones(n, 1);
y0 = -1.6759 * n;
K = 3;
mesh_size = 8;
num_mesh_refine = 7;
max_iter = 35;
verbose = 1;

ddogs = DeltaDOGS;
ddogs.initial(n, bounds, surrogate_, func_eval, xE, y0, K, num_mesh_refine, mesh_size, max_iter, verbose);
ddogs.DeltaDogsOptimize;
