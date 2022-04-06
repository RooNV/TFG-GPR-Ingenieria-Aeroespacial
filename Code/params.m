% Esta función devuelve las funciones necesarias para la implementación de
% diferentes tipos de kernel
% TFG Aero Rocío Navarro Villarino

function [kfn]  = params(tipo)

if isempty(tipo) 
  error('Introduce el tipo de función kernel');
end

    switch tipo
        case 'SE'
            %Kernel cuadrático. (SE)
            kfn = @(x,z,H) (H(1)^2)*exp(-(pdist2(x,z).^2)/(2*H(2)^2));
    %         kfn = @(x,z,H) (H(1)^2)*exp(-sq_dist(x',z')/(2*H(2)^2)); 
        case 'LIN'
            %Kernel Lineal (LIN) 
            kfn = @(x,z,H) (H(1)^2)*kernel_lin(x',z'); 
        case 'PER'
            %Kernel Periódico (PER) 
%             kfn = @(x,z,H) (exp(H(1))^2)*exp(-(2/H(2)^2).*sin(pi/H(3).*kernel_per(x',z')).^2); 
            kfn = @(x,z,H) (H(1)^2)*exp(-(2/H(2)^2)*sin((pi/H(3)).*kernel_per(x',z')).^2); 
        case 'LP'
            %Locally Periodic Kernel (LP)
%             exp para la prediccion 
              kfn = @(x,z,H) (exp(H(1))^2)*exp(-sq_dist(x',z')/(2*H(2)^2)).*exp(-(2/H(3)^2).*sin(pi/H(4).*kernel_per(x',z')).^2);   
%             kfn = @(x,z,H) (H(1)^2)*exp(-sq_dist(x',z')/(2*H(2)^2)).*exp(-(2/H(3)^2).*sin(pi/H(4).*kernel_per(x',z')).^2);
           
        case 'MX1'
            %LIN X SE
            kfn = @(x,z,H) (H(1)^2)*exp(-(pdist2(x,z).^2)/(2*exp(H(2))^2)).*kernel_lin(x',z');
        case 'MX4'
            %SE X SE
%             exp para la prediccion
            kfn = @(x,z,H) (exp(H(1))^2)*(exp(H(3))^2)*exp(-(pdist2(x,z).^2)/(2*H(2)^2)).*exp(-(pdist2(x,z).^2)/(2*H(4)^2));
%             kfn = @(x,z,H) (H(1)^2)*(H(3)^2)*exp(-(pdist2(x,z).^2)/(2*H(2)^2)).*exp(-(pdist2(x,z).^2)/(2*H(4)^2));
    end

end

% case 'POL'
%    %Kernel Polinómico (POL) 
% 	kfn = @(x,z,H) (H(1)^2)*(kernel_lin(x',z') + H(2)).^H(3); 
% 	kfn = @(x,z,H) (H(1)^2)*(kernel_lin(x',z') + 1).^H(2); 
% case 'MX2'
%     kfn = @(x,z,H) (H(1)^2)*kernel_lin(x',z').*(kernel_lin(x',z') + 1).^H(2);
% case 'MX3'
% 	kfn = @(x,z,H) (H(1)^2)*(kernel_lin(x',z') + (kernel_lin(x',z') + 1).^H(3));
    
    
    