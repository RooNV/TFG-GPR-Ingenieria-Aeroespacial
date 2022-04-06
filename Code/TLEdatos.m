function [OE] = TLEdatos(TLE)
%% El objetivo de esta función: 
% Conseguir todos los parametros a partir de la información del TLE.

%% Línea 1:
mu = 398618.0;  % [km3/s2] Parámetro gravitacional estándar de la Tierra.
line1 = TLE(1,:,1);

% The TLE's Epoch indicates the UTC time when its orbit elements were true.
OE.epoch_year = str2double(line1(19:20));     % Epoch year
if (OE.epoch_year < 22)
    OE.epoch_year = OE.epoch_year + 2000;
else
    OE.epoch_year = OE.epoch_year + 1900;
end
OE.epoch_dayOfYear = str2double(line1(21:32));         % Epoch day
OE.epoch = datenum([OE.epoch_year,1,OE.epoch_dayOfYear,0,0,0]);% (datenum) Epoch


%% Línea 2:
line2 = TLE(2,:,1);

% [deg] Orbit Inclination (i) is the angle the satellite's orbit 
% plane makes with the Earth's equatorial plane.  Allowed values are
% [0,180].  Between 0 and 90 degrees, the satellite has prograde motion,
% viewd from a point north of the equator; the satellite moves
% counterclockwise around the Earth.  Between 90 and 180 degrees, the
% satellite has retrograde motion; the satellite moves clockwise
% around the Earth.  
OE.i_deg = str2double(line2(9:16));

% [deg] Right Ascension of Ascending Node is the geocentric Right Ascension
% of the satellite as it intersects the Earth's equatorial plane traveling 
% northward (ascending). Allowed values are [0,360].
OE.Omega_deg = str2double(line2(18:25));    

% [unitless] Eccentricity is the ratio of the orbits focus distance to the
% semi-major axis.  For a circle, the ratio is 0.  As the ellipse becomes
% more elliptical, the eccentricity increases.  Allowed values [0,1).
OE.e = str2double(['0.' line2(27:33)]);     

% [deg] Argument of Perigee (omega or w) is the angle that lies in the
% satellite orbit plane from the Ascending Node (Omega or W) to the perigee
% point (p) along the satellite's direction of travel.
OE.omega_deg = str2double(line2(35:42));    

% [deg] Mean Anomaly parameterizes the location of the satellite on its
% orbit at the time of the TLE epoch. Allowed values [0 360].
%        M(t) = M0 + n(t-t0)
% where
%    M(t) - the Mean Anomaly at time t
%    M0   - a reference Mean Anomaly at time t=0
%    n    - the satellites Mean Motion (orbits/day)
%    t    - time
%    t0   - time of the reference Mean Anomaly
OE.M_deg = str2double(line2(44:51));        

% [km] semi major axis
OE.a_km = ( mu/(str2double(line2(53:63))*2*pi/86400)^2 )^(1/3);  

% [orbits/day] Mean Motion is the number of times the satellite orbits the
% Earth in exactly 24 hours (one solar day).  The theoretical limits are
% between 0 and 17 orbits per day.
%        a = G*Me/(2*pi*n)^2
% where
%        n  - Mean Motion
%        a  - semi-major axis
%        G  - universal gravitational constant
%        Me - mass of the Earth
OE.n_orbits_per_day = line2(53:63);         

end