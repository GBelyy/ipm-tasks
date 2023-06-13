clc
clear all
close all

%% memory allocation
dt = 0.1;                % time step, sec
t = 0 : dt : 500;     % time grid, sec
N = length(t);

params.dt = dt;

%% initial conditions
% physical params
params.earthRadius = 6378.137;     % earth raduis, km
params.earthGM = 3.986004415e5;  % gravitational parameter of the Earth, km^3 / sec^2
params.J2 = 1.08262668e-3;        % J2 coefficient
params.J = [0.52, 0, 0;
            0, 0.58, 0;
            0, 0, 0.71];
params.invJ = inv(params.J);
params.theta = deg2rad(11.7); 
params.torque_k = 3.2; % torque moment
params.RMM = 0*[rand;rand;rand];
params.magMuEarth = 7.812e6; % geomagnetic constant, km^3 кг с^-2 А^-1
params.omegaEarth = 7.29e-5;
params.mass = 30;

% model params
params.gravOption = '2bp';             % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'quat';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'

params.vizualize = true;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 600;        % semi-major axis, km
r0(2) = deg2rad(12);                     % inclination, rad
r0(3) = deg2rad(30);                     % RAAN, rad
r0(4) = 0;                               % eccentricity
r0(5) = deg2rad(0);                      % argument of perigee, rad
r0(6) = deg2rad(0.5);                      % true anomaly, rad
r0 = oe2st(r0, params);

% initial conditions for angular parameters
ang0 = zeros(3,1);
ang0(1) = 0;
ang0(2) = pi/6;
ang0(3) = pi/4;
ang0 = eulerAng2quat(ang0);

% initial conditions for angular velocity
omega0 = zeros(3,1);
omega0(1) = 0.01;
omega0(2) = 0.03;
omega0(3) = 0.04;

x0 = formInitConds(r0, ang0, omega0);
x = zeros(size(ang0,1) + 9, N);
x(:, 1) = x0;

% control parameters
params.controlType = 'no';

%sensors params
params.mag.orientReal  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.mag.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.mag.Expect = 0;
params.mag.sigma = 1e-9;

params.sol.orientReal  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.sol.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.sol.Expect = 0;
params.sol.sigma = deg2rad(0.01);

params.asr.orientReal  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.asr.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.asr.Expect = 0;
params.asr.sigma = 0.001;

params.star.orientReal  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.star.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.star.Expect = 0;
params.star.sigma = 1e-5;

% parameters of Kalman filter
params.filtertype = 'star+ars'; % available settings: 'star+ars', 'ars+mag+sol'

%! after changing filtertype redefine Q, R, P0 matrices
erQ2 = 1e-14;
erW2 = 1e-10;
params.Q = diag([erQ2, erQ2, erQ2, erW2, erW2, erW2]);
sigmaQ2 = 1e-10;
sigmaW2 = 1e-6;

params.R = diag([sigmaQ2, sigmaQ2, sigmaQ2, sigmaW2, sigmaW2, sigmaW2]);
P0 = diag([params.star.sigma^2, params.star.sigma^2, params.star.sigma^2, ...
                 params.asr.sigma^2, params.asr.sigma^2, params.asr.sigma^2]);

%% sensor initialization
sensors = struct();
sensors.mag = Indicator('magnetometer', params.mag);
sensors.sol = Indicator('sunsensor', params.sol);
sensors.ars = Indicator('ars', params.sol);
sensors.star = Indicator('star', params.star);

kalmanData = zeros(7, N-1);
deltaKalman = zeros(7, N-1);
for i = 1:N - 1
    if i == 1
        X_prev = x(:, 1);
        P_prev = P0;
    end
    k1 = rhsAltitude(x(:, i), t(i), params);
    k2 = rhsAltitude(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rhsAltitude(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rhsAltitude(x(:, i) + dt * k3, t(i) + dt, params);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
    
    [X_est, P_est] = calcKalman(t(i), t(i + 1), X_prev, P_prev, params, sensors);
    kalmanData(:, i) = X_est(7:13);
    deltaKalman(:, i) = X_est(7:13) - x(7:13, i+1);
    
    X_prev = x(:,i+1);
    P_prev = P_est;
end

%% visualization
if params.vizualize
    figure
    subplot(2,2,1);
    hold on
    grid on
    xlabel('t, c')
    ylabel('Integrated q')
    plot(t, x(7, :), 'black')
    plot(t, x(8, :), 'red')
    plot(t, x(9, :), 'blue')
    plot(t, x(10, :), 'green')
    legend({'q_0', 'q_1', 'q_2', 'q_3'})
    
    subplot(2,2,2);
    hold on
    grid on
    xlabel('t, c')
    ylabel('Integrated \omega, rad/c')
    plot(t, x(11, :), 'black')
    plot(t, x(12, :), 'red')
    plot(t, x(13, :), 'blue')
    legend({'w_1', 'w_2', 'w_3'})
    
    subplot(2,2,3);
    hold on
    grid on
    xlabel('t, c')
    ylabel('Kalman q')
    plot(t(1:end-1), kalmanData(1, :), 'black')
    plot(t(1:end-1), kalmanData(2, :), 'red')
    plot(t(1:end-1), kalmanData(3, :), 'blue')
    plot(t(1:end-1), kalmanData(4, :), 'green')
    legend({'q_0', 'q_1', 'q_2', 'q_3'})
    
    subplot(2,2,4);
    hold on
    grid on
    xlabel('t, c')
    ylabel('Kalman \omega, rad/c')
    plot(t(1:end-1), kalmanData(5, :), 'black')
    plot(t(1:end-1), kalmanData(6, :), 'red')
    plot(t(1:end-1), kalmanData(7, :), 'blue')
    legend({'w_1', 'w_2', 'w_3'})

    figure
    subplot(2,1,1);
    hold on
    grid on
    xlabel('t, c')
    ylabel('Delta q')
    plot(t(1:end-1), deltaKalman(1, :), 'black')
    plot(t(1:end-1), deltaKalman(2, :), 'red')
    plot(t(1:end-1), deltaKalman(3, :), 'blue')
    plot(t(1:end-1), deltaKalman(4, :), 'green')
    legend({'q_0', 'q_1', 'q_2', 'q_3'})
    
    subplot(2,1,2);
    hold on
    grid on
    xlabel('t, c')
    ylabel('Delta \omega, rad/c')
    plot(t(1:end-1), deltaKalman(5, :), 'black')
    plot(t(1:end-1), deltaKalman(6, :), 'red')
    plot(t(1:end-1), deltaKalman(7, :), 'blue')
    legend({'w_1', 'w_2', 'w_3'})
end

%% auxilary functions
function x0 = formInitConds(r0, ang0, omega0)
    [m, n] = size(ang0);
    x0 = zeros(9+m, n);
    
    x0(1:6, 1) = r0;
    x0(7:6+m, 1:n) = ang0;
    x0(6+m+1:9+m, 1) = omega0;
end