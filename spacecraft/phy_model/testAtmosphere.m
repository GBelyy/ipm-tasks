clc
params.earthRadius = 6378.137;     % earth raduis, km
params.earthGM = 3.986004415e5;    % gravitational parameter of the Earth, km^3 / sec^2
params.magMuEarth = 7.812e6;       % geomagnetic constant, km^3 kg с^-2 А^-1
params.omegaEarth = 7.29e-5;
params.overallSize = 1;            % overall size of spacecraft, m
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 400;        % semi-major axis, km
r0(2) = deg2rad(60);                     % inclination, rad
r0(3) = deg2rad(0);                      % RAAN, rad
r0(4) = 0;                               % eccentricity
r0(5) = deg2rad(0);                      % argument of perigee, rad
r0(6) = deg2rad(0);                      % true anomaly, rad
r0 = oe2st(r0, params);

r = r0(1:3);
v = r0(4:6);
altitude = norm(r) - params.earthRadius;
rho = calcAtmDensity(r, params);

ang0 = zeros(3,1);
ang0(1) = 0;
ang0(2) = 0;
ang0(3) = 0;
ang0 = eulerAng2quat(ang0);

wEarth = [0;0;1] * params.omegaEarth;
vAtm = v - cross(wEarth, r);
ev = vAtm/ norm(vAtm);

vUnit = v / norm(v);
D = quat2dcm(ang0');
%% spacecraft model
normals = [1, -1,  0,  0,  0,  0;
           0,  0,  1, -1,  0,  0;
           0,  0,  0,  0,  1, -1];

x = zeros(1,6);
y = zeros(1,6);
z = zeros(1,6);
for i = 1:6
    n_i = D' * normals(:,i);
    x(i) = n_i(1);
    y(i) = n_i(2);
    z(i) = n_i(3);
end

%hold on
%quiver3(zeros(1,6), zeros(1,6), zeros(1,6) ,x,y,z, 'r')
%quiver3(0, 0, 0 , ev(1), ev(2), ev(3), 'b')


%% visualization
figure 
hold on
quiver3(r(1) * ones(1,6), r(2) * ones(1,6), r(3) * ones(1,6), x, y, z, 'r')
quiver3(r(1), r(2), r(3), ev(1), ev(2), ev(3), 'b')
quiver3(r(1), r(2), r(3), vUnit(1), vUnit(2), vUnit(3), 'g')
grid on
legend({'normals', 'e_v', 'v'})

%%
[F,M] = calcAtmResistance(r,v, D, 0, params);
F
M
