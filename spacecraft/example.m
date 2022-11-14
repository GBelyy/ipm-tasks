clc
clear all
close all

%% memory allocation
dt = 60;                % time step, sec
t = 0 : dt : 3600 * 24;  % time grid, sec
N = length(t);
x = zeros(6, N);

%% initial conditions
% constant coefficients
params.earthRadius = 6378137;     % earth raduis, m
params.earthGM = 3.986004415e14;  % gravitational parameter of the Earth, m^3 / sec^2
params.J2 = 1.08262668e-3;        % J2 coefficient

% model params
params.gravPption = 'j2';             % type of dynamics for calculation, 2bp or j2
params.angParametriztion = 'euler';

% system params
params.checkIntegrals = 0;
params.vizualize = 0;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 600 * 1000; % semi-major axis, m                 
r0(2) = deg2rad(82);                     % inclination, rad                   
r0(3) = deg2rad(270);                    % RAAN, rad           
r0(4) = 0.001;                           % eccentricity                    
r0(5) = deg2rad(0);                      % arguement of perigee, rad          
r0(6) = deg2rad(0);                      % true anomaly, rad
r0 = oe2st(r0, params);

% initial conditions for angular parameters
ang0 = zeros(3,1);
ang0(1) = 0;
ang0(2) = pi/2;
ang0(3) = pi/4;

% initial conditions for angular velocity
omega0 = zeros(3,1);
omega0(1) = 0;
omega0(2) = 0;
omega0(3) = 0.01;

x0 = formInitConds(r0, ang0, omega0);

%% integration
for i = 1:N - 1
    k1 = rightSide(x(:, i), t(i), params);
    k2 = rightSide(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rightSide(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rightSide(x(:, i) + dt * k3, t(i) + dt, params);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
end

%% check first integrals
if params.checkIntegrals
    h = zeros(1, N);
    c = zeros(1, N);
    f = zeros(1, N);
    potEnergy = zeros(1, N);
    
    for i = 1:N
        integrals = calcIntegrals(x(:, i), params);
        h(i) = integrals.h;
        c(i) = norm(integrals.c);
        f(i) = norm(integrals.f);
        
        potEnergy(i) = params.earthGM / norm(x(1:3, 1));
    end
    
    figure
    xlabel('time, sec')
    ylabel('value')
    plot(t, h)
    title('Energy integral')
    grid on
    
    figure
    xlabel('time, sec')
    ylabel('value')
    plot(t, c)
    title('Norm of kinetic moment')
    grid on
    
    figure
    xlabel('time, sec')
    ylabel('value')
    plot(t, f)
    title('Norm of Laplas vector')
    grid on
end

%% vizualization 
if params.vizualize
    figure
    grid on
    plot3(x(1,:), x(2,:), x(3,:))
    xlabel('x, m')
    ylabel('y, m')
    zlabel('z, m')
    grid on
end

%% auxilary functions
function x0 = formInitConds(r0, ang0, omega0)
    n = size(ang0, 2);
    x0 = zeros(12, n);

    x0(1:6, 1) = r0;
    x0(7:9, 1:n) = ang0;
    x0(10:12, 1) = omega0;
end