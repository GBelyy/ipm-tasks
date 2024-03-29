clc
clear all

%% initial conditions
% parameters of physical model
params.m = 0.1;                           % mass of rigid body, kg
params.l = 0.1;                           % half of length of rigid body, m
params.g = 9.86;                          % gravity acceleration, m/s2
params.J = [2, 0, 0;
            0, 2, 0
            0, 0, 10];                     % inertial tensor 
params.invJ = inv(params.J);
params.angParametrization = 'quat';       % type of angular parametrization: 'euler', 'cosineMat', 'quat'
params.integralsEvent     = 'lagrange';    % parameter shows integrals to check: 'euler', 'lagrange'

% operational flags
params.vizualize = 0;
params.checkIntegrals = 1;

% initial conditions for angular parameters
angles = zeros(3, 1);
angles(1) = 0;
angles(2) = pi/4;
angles(3) = 0;

ang0 = eulerAng2quat(angles);
% ang0 = angles;

% A0 = eulerAng2cosineMat(angles);
% ang0 = [A0(:, 1); A0(:, 2); A0(:, 3)];

% initial conditions for angular velocity
omega0 = zeros(3, 1);
omega0(1) = 0;
omega0(2) = 0.1;
omega0(3) = 0;

x0 = formInitConds(ang0, omega0);

%% memory allocation
dt = 1e-6;             % time step, sec
t = 0 : dt : 1;        % time grid, sec
N = length(t);

[m, n] = size(x0);
x = zeros(m, N);
x(1:m, 1) = x0;

%% integration
for i = 1:N - 1
    k1 = rhsAngle(x(:, i), t(i), params);
    k2 = rhsAngle(x(:, i) + dt / 2 * k1, t(i) + dt / 2, params);
    k3 = rhsAngle(x(:, i) + dt / 2 * k2, t(i) + dt / 2, params);
    k4 = rhsAngle(x(:, i) + dt * k3, t(i) + dt, params);
    x(1:m, i + 1) = x(1:m, i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
end

%processing of state vector (quat's normalization, matrix ortogonalization and angles bounding)
%x = processStates(x, params);

%% check integrals
if params.checkIntegrals
    integrals = calcIntegralsAngle(x, params);
    figure
    xlabel('time, sec')
    ylabel('value')
    plot(t, integrals.energy - integrals.energy(1))
%     ylim([min(integrals.energy) - 1 max(integrals.energy) + 1])
    title('Energy')
    grid on

    switch lower(params.integralsEvent)
        case 'euler'
            figure
            xlabel('time, sec')
            ylabel('value')
            plot(t, integrals.K - integrals.K(1))
%             ylim([min(integrals.K)-1 max(integrals.K)+1])
            title('Norm of kinetic moment')
            grid on
        case 'lagrange'     
            figure
            xlabel('time, sec')
            ylabel('value')
            plot(t, integrals.Ke3 - integrals.Ke3(1))
%             ylim([min(integrals.Ke3) - 1 max(integrals.Ke3)+1])
            title('Projection on e3 of kinetic moment')
            grid on
        
            figure
            xlabel('time, sec')
            ylabel('value')
            plot(t, integrals.Kz - integrals.Kz(1))
%             ylim([min(integrals.Kz) - 1 max(integrals.Kz) + 1])
            title('Projection on oz of kinetic moment')
            grid on
    end
end

%% vizualization
if params.vizualize

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
                angles(:, i) = quat2eulerAng(quat(:, i));
            end
        case 'euler'
            for i = 1:N
                angles(:, i) = x(1:3, i);
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
    title(['Quaternion in ', params.angParametrization, ' parametrization'])
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
    title(['Euler angles in ', params.angParametrization, ' parametrization'])
    xlabel('time, sec')
    ylabel('angle,rad')
    legend('\phi','\theta', '\psi')
    grid on
    hold off
    
    k = 4 + int8(isequal(params.angParametrization,'quat')) + 6 * int8(isequal(params.angParametrization,'cosineMat'));
    figure
    grid on
    plot(t, reshape(x(k, :), [1, N]))
    hold on
    plot(t, reshape(x(k + 1, :), [1, N]))
    plot(t, reshape(x(k + 2, :), [1, N]))
    title(['Angular velocity in ', params.angParametrization, ' parametrization'])
    xlabel('time, sec')
    ylabel('omega component, rad/sec')
    legend('p','q', 'e')
    grid on
    hold off

end

%% additional functions
function x0 = formInitConds(ang0, omega0)
    [m, n] = size(ang0);
    x0 = zeros(6, n);

    x0(1:m, 1:n) = ang0;
    x0(m + 1:m + 3, 1) = omega0;
end
function stateProcessed = processStates(x, params)
    % funcion ortogonolize cosine matrices and norm quaternions
    sizes = size(x);
    N = sizes(end);
    stateProcessed = x;
    for i = 1:N
        switch (params.angParametrization)
            case 'cosinemat'
                A = orthogonalize(x(1:3, 1:3, i));
                stateProcessed(1:3, 1:3, i) = A;
            case 'quat'
                quat = x(1:4, 1, i) / norm(x(1:4, 1, i));
                stateProcessed(1:4, i) = quat;
            case 'euler'
                angles = mod(x(1:3, 1, i), 2 * pi);
                stateProcessed(1:3, 1, i) = angles;
        end
    end
end