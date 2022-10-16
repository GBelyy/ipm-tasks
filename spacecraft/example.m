clc
clear all
close all

%% memory allocation
dt = 60;             % time step, sec
t = 0 : dt : 86400; % time grid, sec
N = length(t);
x = zeros(6, N);

%% initial conditions
params.earthRadius = 6378137;     % earth raduis, m
params.earthGM = 3.986004415e14;  % gravitational parameter of the Earth, m^3 / sec^2
params.J2 = 1.08262668e-3;        % J2 coefficient
params.option = 'j2';             % type of dynamics for calculation, 2bp or j2

x0 = zeros(6,1);
x0(1) = params.earthRadius + 600 * 1000; % semi-major axis, m                 
x0(2) = deg2rad(82);                     % inclination, rad                   
x0(3) = deg2rad(270);                    % RAAN, rad           
x0(4) = 0.001;                           % eccentricity                    
x0(5) = deg2rad(0);                      % arguement of perigee, rad          
x0(6) = deg2rad(0);                     % true anomaly, rad

x(:, 1) = oe2st(x0, params);

%% integration
for i = 1:N - 1
    k1 = rightSide(x(:, i), t(i), params);
    k2 = rightSide(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rightSide(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rightSide(x(:, i) + dt * k3, t(i) + dt, params);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
end

%% check first integrals
h = zeros(1, N);
c = zeros(1, N);
f = zeros(1, N);
potEnergy = zeros(1, N);

for i = 1:N
    integrals = calcIntegrals(x(:, i), params);
    h(i) = integrals.h;
    c(i) = norm(integrals.c);
    f(i) = norm(integrals.f);
    
    potEnergy(i) = params.earthGM / norm(x(1:3));
end

figure
hold on
grid on
xlabel('time, sec')
ylabel('value')
plot(t, h)
title('Energy integral')

figure
hold on
grid on
xlabel('time, sec')
ylabel('value')
plot(t, c)
title('Norm of kinetic moment')

figure
hold on
grid on
xlabel('time, sec')
ylabel('value')
plot(t, f)
title('Norm of Laplas vector')

%% vizualization 
figure
grid on
plot3(x(1,:), x(2,:), x(3,:))
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')


%%