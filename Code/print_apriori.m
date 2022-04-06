%% generate data
clear all; close all;
L = 1;
xs = (-10:0.1:10)'; %test data
ns = length(xs);
keps = 1e-7;

%define mean and kernel functions
muFn = @(x) 0*x(:).^2;
% Kernel Cuadrático-exps (SE), Kernel Lineal (LIN), Kernel Periódico (PER), 
% Kernel Polinómico (POL), Kernel mixto (MX).
tipo = 'LP'
Kfn = params(tipo);


%% plot sampled functions from the prior

param = [1 1 1 6];
figure; hold on

model = struct('mu', muFn(xs), 'Sigma',  Kfn(xs, xs, param) );
fs1 = gauss_sample(model, 1);

model = struct('mu', muFn(xs), 'Sigma',  Kfn(xs, xs, param) );
fs2 = gauss_sample(model, 1);

model = struct('mu', muFn(xs), 'Sigma',  Kfn(xs, xs, param));
fs3 = gauss_sample(model, 1);

plot(xs, fs1, xs, fs2, xs, fs3, 'linewidth', 1)