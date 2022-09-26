clc
clear all
close all

%% memory allocation
dt = 0.01;
t = 0:dt:100;
N = length(t);
x = zeros(2, N);

%% initial conditions
x(1, 1) = pi/4;
x(2, 1) = 0.1;

params = struct();
params.g = 9.8;
params.mass = 0.1;
params.l = 0.5;
params.b = 0.1;
params.option = 'control';

params.k_phi = 1;
params.k_omega = 1;

params.x_ref = [pi, 0];

%% integration

for i = 1:N - 1

    k1 = rightSide(x(:,i), t(i), params);
    k2 = rightSide(x(:,i) + dt/2*k1, t(i) + dt/2, params);
    k3 = rightSide(x(:,i) + dt/2*k2, t(i) + dt/2, params);
    k4 = rightSide(x(:,i) + dt*k3, t(i) + dt, params);
    x(:, i + 1) = x(:, i) + dt/6*(k1 + k2*2 + k3*2 + k4);
end

%% check energy balance
energy = zeros(1, N);
kineticEnergy = params.mass * params.l^2 /2 * x(2, :) .^2;
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

if strcmp(params.option, 'control')
    phi_label = 'delta phi, rad';
    omega_label = 'delta omega, rad/s';
    phi = x(1, :) - params.x_ref(1);
    omega = x(2, :) - params.x_ref(2);
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