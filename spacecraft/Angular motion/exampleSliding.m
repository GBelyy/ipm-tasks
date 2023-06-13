clc
clear all
close all

%% memory allocation
dt = 1;                % time step, sec
t = 0 : dt :1600;     % time grid, sec
N = length(t);

params.dt = dt;

%% initial conditions  
params.earthRadius = 6378.137;     % earth raduis, km
params.earthGM = 3.986004415e5;  % gravitational parameter of the Earth, km^3 / sec^2
params.J2 = 1.08262668e-3;        % J2 coefficient
params.J = [2, 0, 0;
            0, 3, 0;
            0, 0, 4];
params.modelJ = spoilTensor(params.J, deg2rad(0.01), 1e-4);
params.invJ = inv(params.J);
params.theta = deg2rad(11.7); 
params.mass = 30;                         % mass of spacecraft, kg
params.torqueMoment = 3.2; % torque moment, A m2

vecRMM = randn(3,1);
params.RMM = vecRMM/norm(vecRMM) * 1e-5; % resudial magnetic moment

params.magMuEarth = 7.812e6; % geomagnetic constant, km^3 kg с^-2 А^-1
params.omegaEarth = 7.29e-5;

% model params
params.gravOption = '2bp';            % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'quat';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'

params.vizualize = true;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 400;        % semi-major axis, km
r0(2) = deg2rad(60);                     % inclination, rad
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
omega0(1) = 0.01;
omega0(2) = 0.05;
omega0(3) = 0.03;
x0 = formInitConds(r0, ang0, omega0);
x = zeros(size(ang0,1) + 9, N);
x(:, 1) = x0;

% control parameters
params.controlType = 'sliding'; % available settings: 'sliding', 'slidingMag'
params.lambda = 1;
params.P = 3 * 1e-3 * eye(3);

%params.lambda0 = 1e-4;
%params.deltaB2 = 0.001;
%params.dotLambda = 0;
%params.lambda = 0.007;
%params.P = 1e-3 * eye(3);
%params.lambda_prev = params.lambda0 * eye(3);

%sensors params
params.mag.orientReal  = [1,0,0;
                          0,1,0;
                          0,0,1];
params.mag.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];

params.sol.orientReal  = [1,0,0;
                          0,1,0;
                          0,0,1];
params.sol.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.sol.Expect = 0;
params.sol.sigma = deg2rad(0.001);
params.mag.Expect = 0;
params.mag.sigma = 1e-6;

discretization = 1;                                              % frequency of control calculation, times per sec
params.controlCalcStep = ceil((1 / discretization) / params.dt); % number of time steps between control calculations
controlDur = 0;

% measurements model
params.measModel = 'triad'; % source of altitude state information: 'ideal', 'triad'

%% sensor initialization
sensors = struct();
sensors.mag = Indicator('magnetometer', params.mag);
sensors.sol = Indicator('sunsensor', params.sol);

%% integration

deltaQuat = zeros(4, N);
deltaOmega = zeros(3, N);
for i = 1:N - 1
    [quat, omega] = measureAltitude(t(i), x(:, i), params, sensors);
    deltaQuat(:, i) = (x(7:10, i) - quat);
    deltaOmega(:, i) = (x(11:13, i) - omega);
    if controlDur == 0
        M_ctrl = calcAngularControl(x(1:3, i), x(4:6, i), quat, omega, t(i), params, sensors);
    end
    controlDur = controlDur + 1;
    if controlDur >= params.controlCalcStep
        controlDur = 0;
    end
    
    k1 = rhsAltitude(x(:, i), t(i), params, M_ctrl);
    k2 = rhsAltitude(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params, M_ctrl);
    k3 = rhsAltitude(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params, M_ctrl);
    k4 = rhsAltitude(x(:, i) + dt * k3, t(i) + dt, params, M_ctrl);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
    x(7:10, i + 1) = x(7:10, i + 1)/norm(x(7:10, i + 1));
end

%% visualization
if params.vizualize
    figure

    subplot(1,2, 1)
    hold on
    grid on
    plot(t, x(11, :), 'blue')
    plot(t, x(12, :), 'green')
    plot(t, x(13, :), 'red')
    legend({'w_{1}', 'w_{2}', 'w_{3}'})
    xlabel('t, c')
    ylabel('w, rad/c')
    grid on

    subplot(1,2,2)
    hold on
    grid on
    plot(t,x(8, :), 'blue')
    plot(t,x(9, :), 'green')
    plot(t,x(10, :), 'red')
    legend({'q_{1}', 'q_{2}', 'q_{3}'})
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