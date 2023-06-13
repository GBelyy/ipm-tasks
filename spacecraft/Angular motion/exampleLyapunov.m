clc
clear all
close all

%% memory allocation
dt = 0.1;             % time step, sec
t = 0 : dt : 800;     % time grid, sec
N = length(t);

params.dt = dt;

%% initial conditions
params.earthRadius = 6378.137;     % earth raduis, km
params.earthGM = 3.986004415e5;    % gravitational parameter of the Earth, km^3 / sec^2
params.J2 = 1.08262668e-3;         % J2 coefficient
params.J = [2, 0, 0;
            0, 3, 0;
            0, 0, 4];
params.modelJ = spoilTensor(params.J, deg2rad(0.01), 1e-4);
params.invJ = inv(params.J);
params.theta = deg2rad(11.7); 
params.mass = 30;                 % mass of spacecraft, kg
params.torqueMoment = 3.2;        % torque moment, A m2

vecRMM = randn(3,1);
params.RMM = vecRMM/norm(vecRMM) * 1e-5; % resudial magnetic moment

params.magMuEarth = 7.812e6; % geomagnetic constant, km^3 кг с^-2 А^-1
params.omegaEarth = 7.29e-5;

% model params
params.gravOption = 'j2';             % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'quat';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'

params.vizualize = true;

% initial conditions for orbital motion
r0 = zeros(6, 1);
r0(1) = params.earthRadius + 400;        % semi-major axis, km
r0(2) = deg2rad(70);                     % inclination, rad
r0(3) = deg2rad(0);                      % RAAN, rad
r0(4) = 0;                               % eccentricity
r0(5) = deg2rad(0);                      % argument of perigee, rad
r0(6) = deg2rad(0);                      % true anomaly, rad
r0 = oe2st(r0, params);

% initial conditions for angular parameters
ang0 = zeros(3,1);
ang0(1) = pi/4;
ang0(2) = pi/4;
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
params.controlType = 'lyapunov';
params.kW = 0.1;
params.kQ = 0.01;
params.omega_prev = [0;0;0];
params.quat_prev = [1;0;0;0];
discretization = 1;                                              % frequency of control calculation, times per sec
params.controlCalcStep = ceil((1 / discretization) / params.dt); % number of time steps between control calculations
controlDur = 0;

%sensors params
params.mag.orientReal  = [1,0,1e-8;
                          1e-8,1,1e-8;
                          0,0,1];
params.mag.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.mag.Expect = 0;
params.mag.sigma = 1e-6;

params.sol.orientReal  = [1,1e-7,0;
                          1e-7,1,1e-7;
                          1e-7,0,1];
params.sol.orientModel  = [1,0,0;
                           0,1,0;
                            0,0,1];
params.sol.Expect = 0;
params.sol.sigma = deg2rad(0.01);

% measurements model
params.measModel = 'triad'; % source of altitude state information: 'ideal', 'triad'


%% sensor initialization
sensors = struct();
sensors.mag = Indicator('magnetometer', params.mag);
sensors.sol = Indicator('sunsensor', params.sol);

%% integration
valueFunc = zeros(1, N);
quatRel = zeros(4, N);
omegaRel = zeros(3, N);

for i = 1:N - 1

    [quat, omega] = measureAltitude(t(i), x(:, i), params, sensors);

    if controlDur == 0
        [M_ctrl, quatRelElem, omegaRelElem, params.omega_prev] = calcAngularControl(x(1:3, i), x(4:6, i), quat, omega, t(i), params, sensors);
    else
        [quatRef, omegaRef] = calcRefMotion(x(1:3,i), x(4:6,i), params);
        quatRelElem = quatmultiply(quatconj(quatRef'), quat')';
        A = quat2dcm(quatRelElem');
        omegaRelElem = omega - A * omegaRef;
        params.omega_prev = omegaRef;
    end

    controlDur = controlDur + 1;
    if controlDur >= params.controlCalcStep
        controlDur = 0;
    end

    valueFunc(i) = 0.5 * dot(omegaRelElem, omegaRelElem) + params.kQ * dot(quatRelElem(2:4), quatRelElem(2:4));
    
    k1 = rhsAltitude(x(:, i), t(i), params, M_ctrl);
    k2 = rhsAltitude(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params, M_ctrl);
    k3 = rhsAltitude(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params, M_ctrl);
    k4 = rhsAltitude(x(:, i) + dt * k3, t(i) + dt, params, M_ctrl);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
    x(7:10, i + 1) = x(7:10, i + 1)/norm(x(7:10, i + 1));
        
    quatRel(:,i) = quatRelElem;
    omegaRel(:, i) = omegaRelElem;
end

%% visualization
if params.vizualize
    
    figure
    hold on
    plot(t(1, 1:end-1), quatRel(2, 1:end-1), 'blue')
    plot(t(1, 1:end-1), quatRel(3, 1:end-1), 'green')
    plot(t(1, 1:end-1), quatRel(4, 1:end-1), 'black')
    legend({'q_{1 rel}', 'q_{2 rel}', 'q_{3 rel}'})
    xlabel('t, c')
    ylabel('q')
    grid on

    figure
    hold on
    plot(t(1, 1:end-1), omegaRel(1, 1:end-1), 'red')
    plot(t(1, 1:end-1), omegaRel(2, 1:end-1), 'blue')
    plot(t(1, 1:end-1), omegaRel(3, 1:end-1), 'green')
    legend({'\delta w_1', '\delta w_2', '\delta w_3'})
    xlabel('t, c')
    ylabel('omega, rad/с')
    grid on

    figure
    xlabel('t, c')
    ylabel('V, rad^2/c^2')
    plot(t(1, 1:end-1), valueFunc(1:end-1), 'red')
    legend({'Value Function'})
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