% TFG Aero Rocío Navarro Villarino
% Función para optimizar los valores de los hiperparámetros
% en este caso con un grid.

function  y_pred = GPR(Xtrain, ftrain, xs, Kfn, W)

%% generate data

    %define mean and kernel functions
    muFn = @(x) 0*x(:).^2;

    %% Hiperparametros
%     Parámetros del grid
%     H_est1 = linspace(1.7781, 1.8812, 20);
%     H_est2 = linspace(115.9196,117.5196, 20);
%     H_est3 = linspace(1.5121, 1.9965, 20);
%     H_est4 = linspace(112.0015,113.1121, 20);
%     Sgn_est = linspace(0.01, 0.01, 1);
%     [rh1, rh2, rh3, rh4, rh5] = ndgrid(H_est1, H_est2, H_est3, H_est4, Sgn_est);
%     
%     x = sym('x');
%     y = sym('y');
%     z = sym('z');
%     t = sym('t');
%     H = [x, y, z, t];
%     
%     SigmaN = @(s) (s^2)*eye(length(Xtrain));
%     K = Kfn(Xtrain, Xtrain, H) + SigmaN;
%     K = matlabFunction(K);
%     result = arrayfun(@(x,y,z,t,s) K(x,y,z,t,s), rh1, rh2, rh3, rh4, rh5, 'UniformOutput',false);
%     
%     fun = zeros(length(H_est1),length(H_est2),length(H_est3),length(H_est4),length(Sgn_est));
%     for i = 1:1:length(H_est1)
%         for j = 1:1:length(H_est2)
%             for p = 1:1:length(H_est3)
%                 for f = 1:1:length(H_est4)
%                     for l = 1:1:length(Sgn_est)
%                         K_subs = cell2mat(result(i,j,p,f,l));
%                         seg = 0.5*log(det(K_subs));
%                         fun(i,j,p,f,l) = - 0.5*ftrain'*inv(K_subs)*ftrain - seg;
%                         if(isinf(fun(i,j,p,f,l)))
%                             fun(i,j,p,f,l) = -Inf(1);
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
%     
% %     Obetención del máximo: La función máximo actúa sobre vectores
%     [~, I] = max(fun(:));
%     [ind1, ind2, ind3, ind4, ind5] = ind2sub(size(fun), I);
%     
% %     Varianles
%     sigma =rh1(ind1,1,1,1,1)
%     L = rh2(1,ind2,1,1,1)
%     sigma2 = rh3(1,1,ind3,1,1)
%     L2 = rh4(1,1,1,ind4,1)
%     sgn = rh5(1,1,1,1,ind5)
%   
%     
%  H = [sigma L sigma2 L2];
    H = [1.8227	116.675 1.7525 113.6965]
    sgn = 0.01;

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
