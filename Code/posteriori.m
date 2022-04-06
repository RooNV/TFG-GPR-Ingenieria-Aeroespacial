% Esta función devuelve las matrices postMU, postCov y Ki para la
% implementación de la a posteriori
% TFG Aero Rocío Navarro Villarino

function [postMu, postCov, Ki]  = posteriori(Kfn, muFn, Xtrain, xs, ftrain, H, sgn)

K = Kfn(Xtrain, Xtrain, H) + sgn^2*eye(length(Xtrain)); % K(x,x) (keps is essential!)
Ks = Kfn(xs, Xtrain, H); %K(x*,x)
Kss = Kfn(xs, xs, H); % K(x*,x*) (keps is essential!)
Ki = inv(K); %O(Ntrain^3)
postMu = muFn(xs) + Ks*Ki*(ftrain - muFn(Xtrain));
postCov = Kss - Ks*Ki*Ks';
end