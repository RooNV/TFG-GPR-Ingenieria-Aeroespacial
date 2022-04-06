% TFG Aero Rocío Navarro Villarino
% Función para optimizar los valores de los hiperparámetros
% en este caso con un grid.

function  y_pred = GPR(Xtrain, ftrain, xs, Kfn, W)

%% generate data

    %define mean and kernel functions
    muFn = @(x) 0*x(:).^2;

    %% Hiperparametros
%     Parámetros del grid
%     H_est1 = linspace(4.6464,4.6464, 1);
%     H_est2 = linspace(0.01212,0.256, 10000);
%     H_est3 = linspace(309.9626,309.9626, 1);
%     Sgn_est = linspace(0.6178, 0.6178, 1);
%     [rh1, rh2, rh3, rh4] = ndgrid(H_est1, H_est2, H_est3, Sgn_est);
%     
%     x = sym('x');
%     y = sym('y');
%     t = sym('t');
%     H = [x, y, t];
%     SigmaN = @(z) (z^2)*eye(length(Xtrain));
%    
%     K = Kfn(Xtrain, Xtrain, H) + SigmaN;
%     
%     K = matlabFunction(K);
%     result = arrayfun(@(x,y,t,z) K(x,y,t,z), rh1, rh2, rh3, rh4, 'UniformOutput',false);
%     
%     fun = zeros(length(H_est1),length(H_est2),length(H_est3),length(Sgn_est));
%     for i = 1:1:length(H_est1)
%         for j = 1:1:length(H_est2)
%             for f = 1:1:length(H_est3)
%                 for p = 1:1:length(Sgn_est)
%                     K_subs = cell2mat(result(i,j,f,p));
%                     K_subs_i = inv(K_subs);
%                     seg_param = 0.5*log(det(K_subs));
%                     pri_param = 0.5*ftrain'*K_subs_i*ftrain;
%                     fun(i,j,f,p) = - pri_param - seg_param;
%                     if(isinf(fun(i,j,f,p)))
%                         fun(i,j,f,p) = -Inf(1);
%                     end
%                 end
%             end
%         end
%     end
%     
%     
% %     Obetención del máximo: La función máximo actúa sobre vectores
%     [~, I] = max(fun(:));
%     [ind1, ind2, ind3, ind4] = ind2sub(size(fun), I);
%     
% %     Varianles
%     sigma =rh1(ind1,1,1,1)
%     L = rh2(1,ind2,1,1)
%     N = rh3(1,1,ind3,1)
%     sgn = rh4(1,1,1,ind4)
%     
%     H = [sigma L N]; 
H = [2.827 0.2446 309.92626]
sgn = 0.6178;

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