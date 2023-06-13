clc
clear all
close all

%% memory allocation
dt = 1;                % time step, sec
t = 0 : dt : 30000;     % time grid, sec
N = length(t);

params.dt = dt;

%% initial conditions
params.earthRadius = 6378.137;           % earth raduis, km
params.earthGM = 3.986004415e5;          % gravitational parameter of the Earth, km^3 / sec^2
params.J2 = 1.08262668e-3;               % J2 coefficient
params.J = [2, 0, 0;
            0, 3, 0;
            0, 0, 4];
params.modelJ = spoilTensor(params.J, deg2rad(0.01), 1e-4);
params.invJ = inv(params.J);
params.theta = deg2rad(11.7); 
params.overallSize = 1;                   % overall size of spacecraft, m
params.torqueMoment = 3.2;                % torque moment, A m2
params.normals = [1, -1,  0,  0,  0,  0;  % unit vector of normals
                  0,  0,  1, -1,  0,  0;  % by spacecraft surfaces in cck
                  0,  0,  0,  0,  1, -1]; 

params.panelNormals = [0.075, -0.075;  % unit vector of normals by solar panels in cck
                       0.075, -0.075; 
                       0.996,  0.996];
params.panelArea = 1;                  % area of solar panels
params.thetaMax = deg2rad(10);         % max angle of solar panel activity
params.mass = 30;                      % mass of spacecraft, kg

vecRMM = randn(3,1);
params.RMM = vecRMM/norm(vecRMM) * 1e-5; % resudial magnetic moment

params.magMuEarth = 7.812e6;             % geomagnetic constant, km^3 кг с^-2 А^-1
params.omegaEarth = 7.29e-5;

% model params
params.gravOption = '2bp';            % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'quat';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'
params.atmOption = true;                 % flag of accounting atmosphere resistance
params.solarOption = true;               % flag of accounting solar pressure
params.vizualize = true;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 600;        % semi-major axis, km
r0(2) = deg2rad(60);                 % inclination, rad
r0(3) = deg2rad(0);                      % RAAN, rad
r0(4) = 0;                               % eccentricity
r0(5) = deg2rad(0);                      % argument of perigee, rad
r0(6) = deg2rad(0);                      % true anomaly, rad
r0 = oe2st(r0, params);

% initial conditions for angular parameters
ang0 = zeros(3,1);
ang0(1) = pi/3;
ang0(2) = pi/6;
ang0(3) = pi/4;
ang0 = eulerAng2quat(ang0);

% initial conditions for angular velocity
omega0 = zeros(3,1);
omega0(1) = 0;
omega0(2) = 0.01;
omega0(3) = 0;
x0 = formInitConds(r0, ang0, omega0);
x = zeros(size(ang0,1) + 9, N);
x(:, 1) = x0;

% control parameters
params.controlType = 'no';

%% sensor initialization
sensors = struct();

%% integration
M_grav = zeros(1,N);
M_atm = zeros(1,N);
M_sol = zeros(1,N);
for i = 1:N - 1
    M_ctrl = calcAngularControl(x(1:3, i), x(4:6, i), x(7:10, i), x(11:13, i), t(i), params, sensors);
    
    k1 = rhsAltitude(x(:, i), t(i), params, M_ctrl);
    k2 = rhsAltitude(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params, M_ctrl);
    k3 = rhsAltitude(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params, M_ctrl);
    k4 = rhsAltitude(x(:, i) + dt * k3, t(i) + dt, params, M_ctrl);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
    x(7:10, i + 1) = x(7:10, i + 1)/norm(x(7:10, i + 1));

    D = quat2dcm(x(7:10, i)');
    [~, M_grav_v, M_atm_v, M_sol_v] = calcExternalMoment(x(:,i), t(i), D, params);

    M_grav(i) = norm(M_grav_v);
    M_atm(i) = norm(M_atm_v);
    M_sol(i) = norm(M_sol_v);
end

%% visualization
if params.vizualize
    quat = x(7:10, :);
    omega = x(11:13, :);

    alt = vecnorm(x(1:3, :)) - params.earthRadius;

    figure
    hold on
    plot(t, quat(1, :), 'red')
    plot(t, quat(2, :), 'blue')
    plot(t, quat(3, :), 'green')
    plot(t, quat(4, :), 'black')
    legend({'q_{1}', 'q_{2}', 'q_{3}'})
    xlabel('t, c')
    ylabel('q')
    grid on

    figure
    hold on
    plot(t, omega(1, :), 'blue')
    plot(t, omega(2, :), 'green')
    plot(t, omega(3, :), 'black')
    legend({'w_1', 'w_2', 'w_3'})
    xlabel('t, rad/c')
    ylabel('omega, 1/с')
    grid on

    figure
    hold on
    xlabel('t, hours')
    ylabel('altitude, km')
    plot(t/3600, alt, 'blue')
    legend('Altitude')
    grid on

    figure
    grid on
    plot3(x(1,:), x(2,:), x(3,:))
    xlabel('x, m')
    ylabel('y, m')
    zlabel('z, m')
    grid on

    figure
    hold on
    plot(t, M_grav, 'red')
    plot(t, M_atm, 'blue')
    plot(t, M_sol, 'green')
    legend({'M_{grav}', 'M_{atm}', 'M_{sol}'})
    xlabel('t, c')
    ylabel('q')
    grid on
end

%% auxilary functions
function x0 = formInitConds(r0, ang0, omega0)
    [m, n] = size(ang0);
    x0 = zeros(9+m, n);
    
    x0(1:6, 1) = r0;
    x0(7:6+m, 1:n) = ang0;
    x0(6+m+1:9+m, 1) = omega0;
end