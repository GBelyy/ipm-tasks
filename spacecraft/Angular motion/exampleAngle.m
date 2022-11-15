clc
clear all
close all

%% memory allocation
dt = 60;                % time step, sec
t = 0 : dt : 3600 * 24;  % time grid, sec
N = length(t);
x = zeros(6, N);

%% initial conditions
params.J = [20, 0, 0;
            0, 20, 0;
            0, 0, 1];
params.angParametriztion = 'euler';   % type of angular parametrization: 'euler', 'cosineMat', 'quat'
params.vizualize = 1;

% initial conditions for angular parameters
ang0 = zeros(3,1);
ang0(1) = pi/6;
ang0(2) = pi/2;
ang0(3) = pi/4;

% initial conditions for angular velocity
omega0 = zeros(3,1);
omega0(1) = 0;
omega0(2) = 0.1;
omega0(3) = 0.1;

x0 = formInitConds(ang0, omega0);
x(:, 1) = x0;

%% integration
for i = 1:N - 1
    k1 = rhsAngle(x(:, i), t(i), params);
    k2 = rhsAngle(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rhsAngle(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rhsAngle(x(:, i) + dt * k3, t(i) + dt, params);
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
    switch lower(params.angParametriztion)
        case 'cosinemat'
            for i = 1:N
                mat_i = x(7:9, 1:3, i);
                quat(1:4, i) = cosineMat2quat(mat_i);
                angles(1:3, i) = cosineMat2eulerAng(mat_i);
            end
        case 'quat'
            for i = 1:N
                quat(:, i) = x(7:10, i);
                angles(:, i) = quat2eulerAng(quat(:, i));
            end
        case 'euler'
            for i = 1:N
                angles(:, i) = x(7:9, i);
                quat(:, i) = eulerAng2quat(angles(:, i));
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
    grid on
    plot(t, angles(1, :))
    hold on
    plot(t, angles(2, :))
    plot(t, angles(3, :))
    xlabel('time, sec')
    ylabel('angle,rad')
    legend('\phi','\theta', '\psi')
    grid on
    hold off

    figure
    grid on
    plot(t, x(10,:))
    hold on
    plot(t, x(11,:))
    plot(t, x(12,:))
    xlabel('time, sec')
    ylabel('omega component, rad/sec')
    legend('p','q', 'e')
    grid on
    hold off
end

%%
function x0 = formInitConds(ang0, omega0)
    n = size(ang0, 2);
    x0 = zeros(6, n);

    x0(1:3, 1:n) = ang0;
    x0(4:6, 1) = omega0;
end