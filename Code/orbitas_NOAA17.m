 %% *Órbita NOAA 17
clear all; close all;

mu = 3.98618e14; % [m3/s2] Earth's geocentric gravitational constant
J2 = 1.08263e-3; % second zone harmonics
R = 6371000; %Radio de la Tierra [m]

%% *Extraer los elementos orbitales del TLE*
% Datos TLE: Hay un total de 3684 líneas
N =3680; %Numero de líneas del fichero a leer. Debe ser par
W = 1; %Intervalo de instantes temporales
TLE = fopen('noaa-17.txt', 'r');

% for i = 1:20
%     xx = fgetl(TLE);
% end;
% N = N -22;

% Se obtienen todas las lineas del TLE
for k = 1:N
    lineas(k,:) = fgetl(TLE);
end

% Se agrupan las líneas de 2 en 2 en TLE_lineas
num_lineas = 0;
for i = 1:2:(N-1)
    num_lineas = num_lineas + 1;
    TLE_filas(:,:,num_lineas) = lineas(i:i+1,:);
end

% Separacion entre instantes temporales
num_puntos = 0;
for i = 1:W:(num_lineas)
    num_puntos = num_puntos + 1;
    TLE_puntos(:,:,num_puntos) = TLE_filas(:,:,i);
end


figure('color','white');
xlabel('ECI x [m]');
ylabel('ECI y [m]');
zlabel('ECI z [m]');
% title('Satellite Orbit in ECI Coordinates');
% pintaTierra();
grid on


for i =  1:num_puntos
   
   % Obtener los datos de las lineas TLE
   TLE = TLE_puntos(:,:,i);
   [OE] = TLEdatos(TLE);

    %% *Matrices de rotación a ECI*
        Rz_Omega(:,:,i) = [ ...
            [cosd(OE.Omega_deg) sind(OE.Omega_deg) 0]; ...
            [-sind(OE.Omega_deg) cosd(OE.Omega_deg) 0]; ...
            [0 0 1]];
        Rx_i(:,:,i) = [ ...
            [1 0 0]; ...
            [0 cosd(OE.i_deg) sind(OE.i_deg)]; ...
            [0 -sind(OE.i_deg) cosd(OE.i_deg)]];
        Rz_omega(:,:,i) = [ ...
            [cosd(OE.omega_deg) sind(OE.omega_deg) 0]; ...
            [-sind(OE.omega_deg) cosd(OE.omega_deg) 0]; ...
            [0 0 1]];

    %% Datos
    Omega_deg(i) = OE.Omega_deg;
    omega_deg(i) = OE.omega_deg;
    i_deg(i) = OE.i_deg; %inclinación
    
    a_m = OE.a_km*1e3; %semi eje mayor de la órbita (m)
    a_mp(i) = a_m;
    e = OE.e; %excentricidad
    e_p(i) = e;
    M_deg = OE.M_deg; %Anomalía media (deg)

    % Anomalía excéntrica
    n_rad_per_s = sqrt(mu/a_m^3);  % [rad/s] mean motion
    n_deg_per_s = rad2deg(n_rad_per_s); % [deg/s] mean motion
    M_rad = deg2rad(M_deg);
    M_rad_p(i)= M_rad;
    E_rad = M_rad; 
    E_rad_p(i) = E_rad;
    dE = 99999;
    eps = 1e-6; % [rad] control precision of Newton's method solution
    while (abs(dE) > eps)
        dE = (E_rad - e * sin(E_rad) - M_rad)/(1 - e * cos(E_rad));
        E_rad = E_rad -  dE;
    end
    
    theta = 2*atan(sqrt((1+e)/(1-e))*tan(E_rad/2));
    p_m = a_m*(cos(theta) - e);
    q_m = a_m*sqrt(1 - e^2)*sin(theta);

    dMdt_rad_per_s = n_rad_per_s;
    dEdt_rad_per_s = dMdt_rad_per_s/(1 - e*cos(E_rad));
    dpdt_m_per_s = -a_m*sin(E_rad)*dEdt_rad_per_s;
    dqdt_m_per_s = a_m*cos(E_rad)*dEdt_rad_per_s*sqrt(1 - e^2);
    E_deg_epoch(i) = rad2deg(E_rad); 

    % Time of epoch 
    Epoch(i) = OE.epoch;
    
 
    %% *Órbitas TLE*
    Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
    Orbit_p = a_m*(cosd(Evals)-e); % [m] orbit positions
    Orbit_q = a_m*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions
    deltaT_s = ((Evals-E_deg_epoch(i)) - e*sind(Evals-E_deg_epoch(i)))/n_deg_per_s; % [s] time since epoch along orbit

    Orbit_ECI = zeros(numel(deltaT_s),3);
    for ipt = 1:size(Orbit_ECI,1)
        r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
        Orbit_ECI(ipt,:) = [inv(Rz_Omega(:,:,i))*inv(Rx_i(:,:,i))*inv(Rz_omega(:,:,i))*r_pq]'; %[Rz_Omega*Rx_i*Rz_omega*r_pq]';
    end

    %% *Posición del punto*
    r_pq = [p_m q_m 0]';
    r_ECI(1,:,i) = [inv(Rz_Omega(:,:,i))*inv(Rx_i(:,:,i))*inv(Rz_omega(:,:,i))*r_pq]';
    
    hold on
%     plot3(r_ECI(:,1,i),r_ECI(:,2,i),r_ECI(:,3,i), '*', 'linewidth', 2)
    
    
end

%% *Órbita del punto 1*
Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
Orbit_p = a_m*(cosd(Evals)-e(1)); % [m] orbit positions
Orbit_q = a_m*sqrt(1 - e_p(1)^2)*sind(Evals); % [m] orbit positions
deltaT_s = ((Evals-E_deg_epoch(1)) - e_p(1)*sind(Evals-E_deg_epoch(1)))/n_deg_per_s; 

Orbit_ECI = zeros(numel(deltaT_s),3);
for ipt = 1:size(Orbit_ECI,1)
    r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
    Orbit_ECI(ipt,:) = [inv(Rz_Omega(:,:,1))*inv(Rx_i(:,:,1))*inv(Rz_omega(:,:,1))*r_pq]'; 
end
% Plot Cartesian Coordinates *Órbita*
% plot3(Orbit_ECI(:,1),Orbit_ECI(:,2),Orbit_ECI(:,3), 'blue', 'linewidth', 2);

%% Punto del Perigeo
p_perigeo = a_m*(cos(0) - e_p(1));
q_perigeo = a_m*sqrt(1 - e_p(1)^2)*sin(0);
r_pq_perigeo = [p_perigeo q_perigeo 0]';

Orbit_ECI_perigeo = [inv(Rz_Omega(:,:,1))*inv(Rx_i(:,:,1))*inv(Rz_omega(:,:,1))*r_pq_perigeo]';
% plot3(Orbit_ECI_perigeo(:,1),Orbit_ECI_perigeo(:,2),Orbit_ECI_perigeo(:,3),'c+', 'linewidth', 1)

%% *Orbital Elements at the inicial time*

T = (2*pi*a_mp(1)^(3/2))/sqrt(mu); % Periodo del primer punto
h = ((T*mu^2/(2*pi))^(1/3))*sqrt(1-e_p(1)^2);
i_rad_1 = deg2rad(i_deg(1));
% e_p(1)
Omega_rad_1 = deg2rad(Omega_deg(1));
omega_rad_1 = deg2rad(omega_deg(1));


%% Comparar la posición de los puntos:
t_p = M_rad_p(1)*T/(2*pi); % segundos desde el perigeo al punto 1

time = zeros(1,(num_puntos-1));
ti = zeros(1,(num_puntos-1));
t_i = zeros(1,(num_puntos-1));
for i=2:1:num_puntos
    % A partir de la anomalía verdadera del primer punto
    % Anomalia media 
    ti(i) = Epoch2seconds(Epoch(i),Epoch(1));
    time(i) = (ti(i) + t_p); %tiempo entre punto2 y perigeo
    M_rad = (2*pi*time(i))/T;

    E_rad = M_rad; % E inicial para la iteración
    dE = 99999;
    eps = 1e-6; % [rad] control precision of Newton's method solution
    while (abs(dE) > eps)
        dE = (E_rad - e_p(i) * sin(E_rad) - M_rad)/(1 - e_p(i) * cos(E_rad));
        E_rad = E_rad - dE;
    end

    theta = 2*atan(sqrt((1+e_p(i))/(1-e_p(i)))*tan(E_rad/2));
    
    p_should = (h^2/(mu*(1+e_p(i)*cos(theta))))*cos(theta);
    q_should = (h^2/(mu*(1+e_p(i)*cos(theta))))*sin(theta);

    r_pq_should = [p_should q_should 0]';
    
    %Regression rate of the ascending node
    derivada_Omega = -1.5*sqrt(mu)*J2*R^2*cos(i_rad_1)/((1-e_p(i)^2)^2*a_mp(i)^(7/2));
    Omega = Omega_rad_1 + derivada_Omega*ti(i);
    
    derivada_omega = (-1.5*sqrt(mu)*J2*R^2/((1-e_p(i)^2)^2*a_mp(i)^(7/2)))*(2.5*sin(i_rad_1)^2-2);
    omega = omega_rad_1 + derivada_omega*ti(i);
    
    Q_Omega(:,:,i) = [ ...
            [cos(Omega) sin(Omega) 0]; ...
            [-sin(Omega) cos(Omega) 0]; ...
            [0 0 1]];
    Q_i(:,:,i) = [ ...
            [1 0 0]; ...
            [0 cos(OE.i_deg) sin(OE.i_deg)]; ...
            [0 -sin(OE.i_deg) cos(OE.i_deg)]];
    Q_omega(:,:,i) = [ ...
            [cos(omega) sin(omega) 0]; ...
            [-sin(omega) cos(omega) 0]; ...
            [0 0 1]];

    Orbit_ECI_should(:,:,i) = [inv(Q_Omega(:,:,i))*inv(Q_i(:,:,i))*inv(Q_omega(:,:,i))*r_pq_should]';
%     plot3(Orbit_ECI_should(:,1,i),Orbit_ECI_should(:,2,i),Orbit_ECI_should(:,3,i),'o')

end; 


%% *Órbita del punto 2*
% for i=2:1:j
%     Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
%     Orbit_p = a_m*(cosd(Evals)-e_p(i)); % [m] orbit positions
%     Orbit_q = a_m*sqrt(1 - e_p(i)^2)*sind(Evals); % [m] orbit positions
%     deltaT_s = ((Evals-E_deg_epoch(i)) - e_p(i)*sind(Evals-E_deg_epoch(i)))/n_deg_per_s; 
% 
%     Orbit_ECI = zeros(numel(deltaT_s),3);
%     for ipt = 1:size(Orbit_ECI,1)
%         r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
%         Orbit_ECI(ipt,:) = [inv(Rz_Omega(:,:,i))*inv(Rx_i(:,:,i))*inv(Rz_omega(:,:,i))*r_pq]';
%     end
% % Plot Cartesian Coordinates *Órbita*
% plot3(Orbit_ECI(:,1),Orbit_ECI(:,2),Orbit_ECI(:,3), 'red');
% end;


%% Error de posición:
error_ejex = zeros(1,(num_puntos));
error_ejey = zeros(1,(num_puntos));
error_ejez = zeros(1,(num_puntos));
for i=2:1:num_puntos
    pto_x1 = Orbit_ECI_should(:,:,i); % Estas son las coordenadas del punto 2, según el punto 1.
    pto_x = r_ECI(:,:,i); %Coordenadas del punto 2 según el TLE.
    
    error_ejex(i) = (pto_x(:,1)-pto_x1(:,1))/1e3;
    error_ejey(i) = (pto_x(:,2)-pto_x1(:,2))/1e3;
    error_ejez(i) = (pto_x(:,3)-pto_x1(:,3))/1e3;
    error_posicion(:,:,i) = [error_ejex(i) error_ejey(i) error_ejez(i)];
    
end;

% Separacion entre instantes temporales
num = 0;
W = 10;
for i = 1:W:(num_puntos)
    num = num + 1;
    errX(num) = error_ejex(i);
    errY(num) = error_ejey(i);
    errZ(num) = error_ejez(i);
end

% PREDICCIÓN/ESTIMACIÓN
% F = 1800/8 - 25;
F = 1600/W - 20; % F = Numero de instantes de entrenamiento: 140 puntos de entrenamiento y 20 de prediccion
tiempo = ti/(3600*24);
Ttrain_gpr = tiempo(1,1:W:F*W)';% Vector temporal de entrenamiento [días]

err = errX;
ftrain_gpr = err(1,1:F)'; % Vector de datos de entrenamiento

% Normalizamos los datos
media = mean(ftrain_gpr);
sigma = std(ftrain_gpr);
ftrain_gpr = (ftrain_gpr - media)/sigma;

Kfn = params('PER'); % Función kernel 
ttest = tiempo(1,1:1:1600)'; % Puntos a estimar

% Predicción Matlab
% gprMdl = Mat_gpr(Ttrain_gpr, ftrain_gpr, Kfn); 
% ypred = predict(gprMdl, ttest);
% 
% figure()
% hold on
% ypred = ypred*sigma + media;
% plot(ttest,ypred,'g','LineWidth',1.5);
% 
% plot(tiempo(1,1401),err(1,F+1),'*k','LineWidth',1.5);
% 
% error_real = error_ejex(1,1:1600);
% plot(ttest,error_real,'k');
% hold off

    
% Prueba con la función GPR creada en este proyecto: el número depende de
% los hyperparametros a optimizar
%  ypred_Ro = gp_reg_grid_3(Ttrain_gpr, ftrain_gpr, ttest, Kfn, W);
ypred_Ro = gp_reg_grid_2(Ttrain_gpr, ftrain_gpr, ttest, Kfn, W);
% ypred_Ro = gp_reg_grid(Ttrain_gpr, ftrain_gpr, ttest, Kfn, W);
% ypred_Ro = gp_reg_grid_0(Ttrain_gpr, ftrain_gpr, ttest, Kfn, W);


% Plot: 
figure()
hold on;

ypred_Ro = ypred_Ro*sigma + media;
plot(ttest,ypred_Ro','Color', [0.3010 0.7450 0.9330],'LineWidth',1.5);

error_real = error_ejex(1,1:1600);
plot(ttest,error_real,'k');

plot(tiempo(1,1401),err(1,F+1),'*k','LineWidth',1.5);

xlabel('Tiempo (días)');
ylabel('Error (km)');
hold off


%% Gráficas del error
% error_ejex = zeros(1,(num_puntos-1));
% error_ejey = zeros(1,(num_puntos-1));
% error_ejez = zeros(1,(num_puntos-1));
% j = 5;
% for i=2:1:num_puntos
%     pto_x1 = Orbit_ECI_should(:,:,i); % Estas son las coordenadas del punto 2, según el punto 1.
%     pto_x = r_ECI(:,:,i); %Coordenadas del punto 2 según el TLE.
%     
%     error_ejex(j) = (pto_x(:,1)-pto_x1(:,1))/1e3;
%     error_ejey(j) = (pto_x(:,2)-pto_x1(:,2))/1e3;
%     error_ejez(j) = (pto_x(:,3)-pto_x1(:,3))/1e3;
%     error_posicion(:,:,i) = [error_ejex(i) error_ejey(i) error_ejez(i)];
%     j = j + 5;
% end;
% 
% figure()
% % subplot(3,1,1);
% bar(error_ejez)
% grid on
% % title('Error de posición eje x')
% xlabel('Muestras')
% ylabel('Error (km)')
% axis auto


% % subplot(3,1,2);
% bar(error_ejey)
% grid on
% % title('Error de posición eje y')
% xlabel('Muestras')
% ylabel('Error (km)')
% axis auto
%     
% % subplot(3,1,3);
% xlim([2 j]);
% bar(error_ejez)
% grid on
% % title('Error de posición eje z')
% xlabel('Muestras')
% ylabel('Error (km)')
% axis auto
    


