clc
clear all
close all

%% memory allocation
dt = 1;                % time step, sec
t = 0 : dt : 3000;     % time grid, sec
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
params.invJ = inv(params.J);
params.mass = 42;
params.theta = deg2rad(11.7); 
params.torque_k = 3.2; % torque moment
params.RMM = 0.5*[rand;rand;rand];
params.magMuEarth = 7.812e6; % geomagnetic constant, km^3 кг с^-2 А^-1
params.omegaEarth = 7.29e-5;

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
omega0(1) = 0.0003;
omega0(2) = 0.0001;
omega0(3) = 0.0004;

x0 = formInitConds(r0, ang0, omega0);
x = zeros(size(ang0,1) + 9, N);
x(:, 1) = x0;

% control parameters
params.controlType = 'no';

%sensors params
params.mag.orientReal  = eye(3);
params.mag.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.mag.Expect = 0;
params.mag.sigma = 1e-6;

params.sol.orientReal  = eye(3);%[1,1e-7,0;
                          %1e-7,1,1e-7;
                          %1e-7,0,1];
params.sol.orientModel  = [1,0,0;
                           0,1,0;
                           0,0,1];
params.sol.Expect = 0;
params.sol.sigma = deg2rad(0.01);

%% sensor initialization
sensors = struct();
sensors.mag = Indicator('magnetometer', params.mag);
sensors.sol = Indicator('sunsensor', params.sol);

triadData = zeros(4, N);
deltaTriad = zeros(4, N);
deltaS = zeros(3, N);
deltaB = zeros(3,N);
cosBS = zeros(1,N);
for i = 1:N - 1
    
    k1 = rhsAltitude(x(:, i), t(i), params);
    k2 = rhsAltitude(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rhsAltitude(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rhsAltitude(x(:, i) + dt * k3, t(i) + dt, params);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
    triadData(:, i) = calcTRIAD(t(i), x(:, i), params, sensors.mag, sensors.sol);
    deltaTriad(:, i) = abs(triadData(:, i) - x(7:10, i));
    
    D = quat2dcm(x(7:10, i)');
    
    JD = 2459580.5 + t(i)/86400;
    sunVec = D * sun(JD)'/norm(sun(JD));
    s = sensors.sol.measure(t(i), x(:,i), params);
    s = s / norm(s);
    deltaS(:,i) = sunVec - s;

    B_i = calcMagneticField(x(1:6,i), t(i), params);
    B_cck = D * B_i;
    B_cck = B_cck/norm(B_cck);
    B = sensors.mag.measure(t(i), x(:,i), params);
    B = B/ norm(B);
    deltaB(:,i) = B_cck - B;

    cosBS(i) = dot(B_cck, sunVec);
end

%% visualization
if params.vizualize
    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('TRIAD')
    plot(t, triadData(1, :), 'black')
    plot(t, triadData(2, :), 'red')
    plot(t, triadData(3, :), 'blue')
    plot(t, triadData(4, :), 'green')
    legend({'q_0', 'q_1', 'q_2', 'q_3'})

    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('q')
    plot(t, x(7, :), 'black')
    plot(t, x(8, :), 'red')
    plot(t, x(9, :), 'blue')
    plot(t, x(10, :), 'green')
    legend({'q_0', 'q_1', 'q_2', 'q_3'})

    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('\Delta')
    plot(t, deltaTriad(1, :), 'black')
    plot(t, deltaTriad(2, :), 'red')
    plot(t, deltaTriad(3, :), 'blue')
    plot(t, deltaTriad(4, :), 'green')
    legend({'q_0', 'q_1', 'q_2', 'q_3'})
    
    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('\Delta s')
    plot(t, deltaS(1, :), 'red')
    plot(t, deltaS(2, :), 'blue')
    plot(t, deltaS(3, :), 'green')
    legend({'s_1', 's_2', 's_3'})

    figure
    hold on
    grid on
    xlabel('t, c')
    ylabel('\Delta B')
    plot(t, deltaB(1, :), 'red')
    plot(t, deltaB(2, :), 'blue')
    plot(t, deltaB(3, :), 'green')
    legend({'B_1', 'B_2', 'B_3'})

    figure
    grid on
    xlabel('t, c')
    ylabel('Angle s')
    plot(t, cosBS, 'red')
    legend({'cos Bs'})
end

%% auxilary functions
function x0 = formInitConds(r0, ang0, omega0)
    [m, n] = size(ang0);
    x0 = zeros(9+m, n);
    
    x0(1:6, 1) = r0;
    x0(7:6+m, 1:n) = ang0;
    x0(6+m+1:9+m, 1) = omega0;
end