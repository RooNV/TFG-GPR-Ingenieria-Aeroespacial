% TFG Aero Rocío Navarro Villarino
% Función para optimizar los valores de los hiperparámetros
% en este caso con un grid.


function  y_pred = GPR(Xtrain, ftrain, xs, Kfn, W)

%% generate data

    %define mean and kernel functions
    muFn = @(x) 0*x(:).^2;

    %% Hiperparametros
%     Parámetros del grid
    H_est1 = linspace(0,2, 50);
    Sgn_est = linspace(0.1, 0.99, 20);
    [rh1, rh3] = ndgrid(H_est1, Sgn_est);
    
    x = sym('x');
    H = [x];
    SigmaN = @(z) (z^2)*eye(length(Xtrain));
   
    K = Kfn(Xtrain, Xtrain, H) + SigmaN;
    
    K = matlabFunction(K);
    result = arrayfun(@(x,z) K(x,z), rh1, rh3, 'UniformOutput',false);
    
    fun = zeros(length(H_est1),length(Sgn_est));
    for i = 1:1:length(H_est1)
        for p = 1:1:length(Sgn_est)
            K_subs = cell2mat(result(i,p));
            K_subs_i = inv(K_subs);
            seg_param = 0.5*log(det(K_subs));
            pri_param = 0.5*ftrain'*K_subs_i*ftrain;
            fun(i,p) = - pri_param - seg_param;
            if(isinf(fun(i,p)))
                fun(i,p) = -Inf(1);
            end
        end
    end
    
    
%     Obetención del máximo: La función máximo actúa sobre vectores
    [~, I] = max(fun(:));
    [ind1, ind3] = ind2sub(size(fun), I);
    
%     Varianles
    sigma =rh1(ind1,1)
    sgn = rh3(1,ind3)
    
    H = [sigma];
%     H = 2.7660*10^(-4);
%     sgn = 0.99;
    [postMu_H,postCov_H,Ki_H] = posteriori(Kfn, muFn, Xtrain, xs, ftrain, H, sgn);


    %% Predicción 
    model_H = struct('mu', postMu_H, 'Sigma', postCov_H);
    y_pred = gauss_sample(model_H, 1);
   
    
end

    function S = gauss_sample(model, n)
        % Returns n samples from a multivariate Gaussian distribution
        % S = AZ + mu

        mu = model.mu;
        Sigma = model.Sigma;

        Sigma = nearestSPD(Sigma);
        A = chol(Sigma, 'lower');
        Z = randn(length(mu), n);
        S = bsxfun(@plus, mu(:), A*Z)';
    end