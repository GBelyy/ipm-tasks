clc
clear all
close all

%% memory allocation
dt = 1;                     % time step, sec
t = 0 : dt : 15000;        % time grid, sec
N = length(t);

params.dt = dt;

%% initial conditions
% physical params
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
params.torqueMoment = 30; % torque moment, A m2

vecRMM = randn(3,1);
params.RMM = vecRMM/norm(vecRMM) * 0.3; % resudial magnetic moment

params.magMuEarth = 7.812e6; % geomagnetic constant, km^3 кг с^-2 А^-1
params.omegaEarth = 7.29e-5;

%sensors params
params.mag.orientReal  = eye(3);%[1,0,1e-8;
                          %1e-8,1,1e-8;
                          %0,0,1];
params.mag.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.mag.Expect = 0;
params.mag.sigma = 1e-8;

% model params
params.gravOption = '2bp';             % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'quat';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'

% control parameters
params.controlType = 'bdot';
params.bdot_coef = 1e7;

params.vizualize = true;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 300;        % semi-major axis, km
r0(2) = deg2rad(10);                     % inclination, rad
r0(3) = deg2rad(30);                     % RAAN, rad
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
omega0(1) = 0.002;
omega0(2) = 0.001;
omega0(3) = 0.004;
%omega0 = randn(3, 1)*1e-2;

x0 = formInitConds(r0, ang0, omega0);
x = zeros(size(ang0,1) + 9, N);
x(:, 1) = x0;

discretization = 1;                                      % frequency of control calculation, times per sec
params.controlCalcStep = ceil((1 / discretization) / params.dt); % number of time steps between control calculations

%% sensor initialization
sensors = struct();
sensors.mag = Indicator('magnetometer', params.mag);

%% integration
controlDur = 0;
for i = 1:N - 1
    if controlDur == 0
        M_ctrl = calcAngularControl(x(1:3, i), x(4:6, i), x(7:10, i), x(11:13, i), t(i), params, sensors);
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
end

%% visualization
if params.vizualize
    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('omega, rad/с')
    plot(t, x(11, :), 'red')
    plot(t, x(12, :), 'blue')
    plot(t, x(13, :), 'green')
    legend({'w_1', 'w_2', 'w_3'})

    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('q')
    plot(t, x(8, :), 'red')
    plot(t, x(9, :), 'blue')
    plot(t, x(10, :), 'green')
    legend({'q_1', 'q_2', 'q_3'})
    
    moduleOmega = zeros(1, N);
    refOmega = zeros(1, N); 
    T = 2 * pi *sqrt((params.earthRadius + 300)^3/params.earthGM);
    for i = 1:N
        moduleOmega(i) = norm(x(11:13, i));
        refOmega(i) = 2 * pi/T * 1.8;
    end

    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('w, rad/с')
    plot(t, moduleOmega, 'red')
    plot(t, refOmega, 'blue')
    legend({'omega', '1.8 omega0'})
end

%% auxilary functions
function x0 = formInitConds(r0, ang0, omega0)
    [m, n] = size(ang0);
    x0 = zeros(9+m, n);
    
    x0(1:6, 1) = r0;
    x0(7:6+m, 1:n) = ang0;
    x0(6+m+1:9+m, 1) = omega0;
end