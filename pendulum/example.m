clc
clear all
close all

%% memory allocation
dt = 0.01;        % time step, sec
t = 0 : dt : 300; % time grid, sec
N = length(t);
x = zeros(2, N);

%% initial conditions
x(1, 1) = pi / 4;                 % initial angle
x(2, 1) = 0.1;                    % initial angular velocity

params = struct();

params.g = 9.8;                   % gravitational acceleration, m/sec^2
params.mass = 0.1;                % mass of pendulum, kg
params.l = 0.5;                   % length of pendulum, m
params.b = 0.1;                   % dumping coefficient, sec^-1

params.option = 'control';        % type of equation's right side

params.koefPhi = 1;               % proportional coefficient of control
params.koefOmega = 1;             % derivative coefficient of control

params.xRef = [pi, 0];            % refrence position

params.stepsForDetermination = 5; % steps for position determination (set 1 to disable discretization)
params.stepsForControl = 5;       % steps for control application (set 1 to disable discretization)

params.sigmaPhi = 0.01;              % angular error (set 0 to disable inaccuracy)
params.sigmaOmega = 0.01;            % angular velocity error (set 0 to disable inaccuracy)

%% integration

sigma_matrix = diag([params.sigmaPhi, params.sigmaOmega]);
x_meas = [0; 0];
cntrl = 0;

for i = 1:N - 1

    if rem(i, params.stepsForDetermination) == 0
        x_meas = x(:, i) + sigma_matrix * randn(2, 1);
    end
    if rem(i, params.stepsForControl) == 0
        cntrl = getControl(x_meas, params);
    end


    k1 = rightSide(x(:, i), t(i), params, cntrl);
    k2 = rightSide(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params, cntrl);
    k3 = rightSide(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params, cntrl);
    k4 = rightSide(x(:, i) + dt * k3, t(i) + dt, params, cntrl);
    x(:, i + 1) = x(:, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
end

%% check energy balance

energy = zeros(1, N);
kineticEnergy = params.mass * params.l^2 /2 * x(2, :) .^ 2;
potentialEnergy = -params.mass * params.g * params.l * cos(x(1, :));

figure
hold on
grid on
xlabel('time, sec')
ylabel('energy, J')
title('Energy balance')
plot(t, kineticEnergy)
plot(t, potentialEnergy)
plot(t, kineticEnergy + potentialEnergy)
legend('Kinetic energy', 'Potential energy', 'Total energy')

%% vizualization
% if option == 'control' plots of phi and omega deviations from reference position will be shown,
% otherwise true phi and omega positions will be plotted.

if strcmp(params.option, 'control')
    phi_label = 'delta phi, rad';
    omega_label = 'delta omega, rad/s';
    phi = x(1, :) - params.xRef(1);
    omega = x(2, :) - params.xRef(2);
else
    phi_label = 'phi, rad';
    omega_label = 'omega, rad/s';
    phi = x(1, :);
    omega = x(2, :);
end

% plot angular relation
figure
hold on
grid on
xlabel('time, sec')
ylabel(phi_label)
plot(t, phi)

% plot angular velocity relation
figure
hold on
grid on
xlabel('time, sec')
ylabel(omega_label)
plot(t, omega)

% plot phase trajectory
figure
hold on
grid on
xlabel('phi, rad')
ylabel('angular velocity, rad/s')
title('Phase trajectory', 'interpreter','latex')
plot(x(1, :), x(2, :))