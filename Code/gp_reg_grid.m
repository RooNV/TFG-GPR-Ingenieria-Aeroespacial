% TFG Aero Rocío Navarro Villarino
% Función para optimizar los valores de los hiperparámetros
% en este caso con un grid.


function y_pred = GPR(Xtrain, ftrain, xs, Kfn, W)

%% generate data

    %define mean and kernel functions
    muFn = @(x) 0*x(:).^2;

    %% Hiperparametros
%     Parámetros del grid
%     H_est1 = linspace(0.0015,0.11111, 30);
%     H_est2 = linspace(2.021,3.2109, 30);
%     Sgn_est = linspace(0.009, 0.01, 2);
%     [rh1, rh2, rh3] = ndgrid(H_est1, H_est2, Sgn_est);
%     
%     x = sym('x');
%     y = sym('y');
%     H = [x, y];
%     SigmaN = @(z) (z^2)*eye(length(Xtrain));
%    
%     K = Kfn(Xtrain, Xtrain, H) + SigmaN;
%     
%     K = matlabFunction(K);
%     result = arrayfun(@(x,y,z) K(x,y,z), rh1, rh2, rh3, 'UniformOutput',false);
%     
%     fun = zeros(length(H_est1),length(H_est2),length(Sgn_est));
%     for i = 1:1:length(H_est1)
%         for j = 1:1:length(H_est2)
%             for p = 1:1:length(Sgn_est)
%                 K_subs = cell2mat(result(i,j,p));
%                 K_subs_i = inv(K_subs);
%                 seg_param = 0.5*log(det(K_subs));
%                 pri_param = 0.5*ftrain'*K_subs_i*ftrain;
%                 fun(i,j,p) = - pri_param - seg_param;
%                 if(isinf(fun(i,j,p)))
%                     fun(i,j,p) = -Inf(1);
%                 end
%             end
%         end
%     end
%     
%     
% %     Obetención del máximo: La función máximo actúa sobre vectores
%     [~, I] = max(fun(:));
%     [ind1, ind2, ind3] = ind2sub(size(fun), I);
%     
% %     Varianles
%     sigma =rh1(ind1,1,1)
%     L = rh2(1,ind2,1)
%     sgn = rh3(1,1,ind3)
%     
%     H = [sigma L];
H = [2.4916	10.2164]
sgn = 0.007;

    try
        [postMu_H,postCov_H,Ki_H] = posteriori(Kfn, muFn, Xtrain, xs, ftrain, H, sgn); 
        model_H = struct('mu', postMu_H, 'Sigma', postCov_H);
        y_pred = gauss_sample(model_H, 1);
    catch
        [postMu_H,postCov_H,Ki_H] = posteriori_ayuda(Kfn, muFn, Xtrain, xs, ftrain, H, sgn); 
        model_H = struct('mu', postMu_H, 'Sigma', postCov_H);
        y_pred = gauss_sample(model_H, 1);
        text = 'Con ayuda'
    end
   
    
end