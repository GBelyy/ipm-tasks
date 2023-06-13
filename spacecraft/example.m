clc
clear all
close all

%% memory allocation
dt = 1e-4;                % time step, sec
t = 0 : dt : 10;  % time grid, sec
N = length(t);

%% initial conditions
% physical params
params.earthRadius = 6378137;     % earth raduis, m
params.earthGM = 3.986004415e14;  % gravitational parameter of the Earth, m^3 / sec^2
params.J2 = 1.08262668e-3;        % J2 coefficient
params.J = [200, 0, 0;
            0, 300, 0;
            0, 0, 400];
params.invJ = inv(params.J);
params.mm = [1;2;3];
params.mu_B = 7.812*10^6;
params.theta = 12*pi\180;
params.omega = 2*pi/(24*3600);

% model params
params.gravOption = 'j2';             % type of dynamics for calculation: '2bp' or 'j2'
params.angParametrization = 'euler';   % type of angular parametrization: 'euler', 'cosinemat', 'quat'

% system params
params.controlType = NaN; 
params.checkIntegrals = false;
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

%% vizualization
if params.vizualize
    figure
    grid on
    plot3(x(1,:), x(2,:), x(3,:))
    xlabel('x, m')
    ylabel('y, m')
    zlabel('z, m')
    grid on

    % calculation of euler angles and quaternions
    quat = zeros(4, N);
    angles = zeros(3, N);
    switch lower(params.angParametrization)
        case 'cosinemat'
            for i = 1:N
                mat_i = [x(1:3, i), x(4:6, i), x(7:9, i)];
                quat(1:4, i) = dcm2quat(mat_i);
                angles(1:3, i) = cosineMat2eulerAng(mat_i);
            end

        case 'quat'
            for i = 1:N
                quat(:, i) = x(1:4, i);
                angles(:, i) = mod(quat2eulerAng(quat(:, i)), 2 * pi);
            end
        case 'euler'
            for i = 1:N
                angles(:, i) = mod(x(1:3, i), 2*pi);
                quat(:, i) = eulerAng2quat(x(1:3, i));
            end
    end

    figure
    grid on
    plot(t, quat(1, :))
    hold on
    plot(t, quat(2, :))
    plot(t, quat(3, :))
    plot(t, quat(4, :))
    xlabel('time, sec')
    ylabel('quat. component')
    legend('q1','q2', 'q3', 'q4')
    grid on
    hold off

    figure
    plot(t, angles(2, :))
    xlabel('time, sec')
    ylabel('angle,rad')
    legend('\theta')
    grid on

    figure
    plot(t, angles(3, :))
    xlabel('time, sec')
    ylabel('angle,rad')
    legend('\psi')
    grid on

    figure
    plot(t, angles(1, :))
    legend('\phi')
    grid on
    hold off

    k = 10 + int8(isequal(params.angParametrization,'quat')) + 6 * int8(isequal(params.angParametrization,'cosineMat'));
    figure
    grid on
    plot(t, x(k, :))
    hold on
    plot(t, x(k + 1, :))
    plot(t, x(k + 2, :))
    title('Angular velocity')
    xlabel('time, sec')
    ylabel('omega component, rad/sec')
    legend('p','q', 'e')
    grid on
    hold off
end

%% calc Jacob's integral
if params.checkIntegrals
    k = 10 + int8(isequal(params.angParametrization,'quat')) + 6 * int8(isequal(params.angParametrization,'cosineMat'));
    jacIntegral = zeros(1, N);
    for i = 1:N
        A = quat2dcm(quat(:, i)');
        E2 = A * [0;1;0];
        E3 = A * [0;0;1];

        omega = x(k:k+2,1);
        r = x(1:3,1);
        jacIntegral(i) = 0.5 * dot(omega, params.J * omega) ...
                         + 3/2 * params.earthGM / norm(r)^3 * dot(E3, params.J *E3) ...
                         - dot(sqrt(params.earthGM / norm(r)) * E2, E2);
    end

    figure
    plot(t, jacIntegral)
    ylim([min(jacIntegral)-1 max(jacIntegral)+1])
    title('Jacobi integral')
    xlabel('time, sec')
    ylabel('integral')
    grid on
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