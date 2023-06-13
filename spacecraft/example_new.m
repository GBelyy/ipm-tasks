clc
clear all
close all

%% memory allocation
dt = 1e-2;                % time step, sec
t = 0 : dt : 1000;        % time grid, sec
N = length(t);

%% initial conditions
% physical params
params.earthRadius = 6378137;     % earth raduis, m
params.earthGM = 3.986004415e14;  % gravitational parameter of the Earth, m^3 / sec^2
params.J2 = 1.08262668e-3;        % J2 coefficient
params.J = [2, 0, 0;
            0, 3, 0;
            0, 0, 4];
params.invJ = inv(params.J);
params.mu_B = 7.812*10^6;

%sensors params
orientModel = eye(3);
orientReal  = [1,0,1e-6;
               1e-6,1,1e-6;
               1e-6,0,1];
% model params
params.gravOption = '2bp';             % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'quat';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'
params.calcControl = 'Bdot';
params.vizualize = 1;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 600 * 1000; % semi-major axis, m
r0(2) = deg2rad(82);                     % inclination, rad
r0(3) = deg2rad(270);                    % RAAN, rad
r0(4) = 0;                               % eccentricity
r0(5) = deg2rad(0);                      % argument of perigee, rad
r0(6) = deg2rad(0);                      % true anomaly, rad
r0 = oe2st(r0, params);

% initial conditions for angular parameters
ang0 = zeros(3,1);
ang0(1) = 0;
ang0(2) = pi/6;
ang0(3) = pi/4;
ang0 = eulerAng2quat(ang0);

% initial conditions for angular velocity
omega0 = zeros(3,1);
omega0(1) = 0;
omega0(2) = 0.01;
omega0(3) = 0.04;

x0 = formInitConds(r0, ang0, omega0);
x = zeros(size(ang0,1) + 9, N);
x(:, 1) = x0;

%% integration
for i = 1:N - 1
    k1 = rightSide(x(:, i), t(i), params);
    k2 = rightSide(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rightSide(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rightSide(x(:, i) + dt * k3, t(i) + dt, params);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
end

%% auxilary functions
function x0 = formInitConds(r0, ang0, omega0)
    [m, n] = size(ang0);
    x0 = zeros(9+m, n);
    
    x0(1:6, 1) = r0;
    x0(7:6+m, 1:n) = ang0;
    x0(6+m+1:9+m, 1) = omega0;
end