function [dot_x] = rightSide(x, t, params)
    [m, n] = size(x);
    dot_x = zeros(m, n);
    M = 0;

    % Вычисление орбитального движения
    switch lower(params.gravOption)
        case '2bp'
            r = x(1:3, 1);
            v = x(4:6, 1);
            dot_x(1:3, 1) = v;
            dot_x(4:6, 1) = - params.earthGM / norm(r)^3 * r;
        case 'j2'
            r = x(1:3, 1);
            v = x(4:6, 1);
            coef = 3/2 * params.J2 * params.earthGM * params.earthRadius^2;
            accelerationJ2 = coef * r / norm(r)^5 * (5 * r(3)^2 / norm(r)^2 - 1) - 2 * coef / norm(r)^5 * [0; 0; r(3)];
            
            dot_x(1:3, 1) = v;
            dot_x(4:6, 1) = - params.earthGM / norm(r)^3 * r + accelerationJ2;
    end

    % Вычисление углового движения
    switch lower(params.angParametriztion)
        case 'euler'
            phi = x(7, 1);
            theta = x(8, 1);
            psi = x(9, 1);
            omega = x(10:12);

            dot_angles(1) = omega(1) * cos(phi) - omega(2) * sin(phi);
            dot_angles(2) = (omega(1) * sin(phi) + omega(2) * cos(phi)) / sin(theta);
            dot_angles(3) = omega(3) - cot(theta) * (omega(1) * sin(phi) + omega(2) * cos(phi));

            dot_omega = inv(params.J) * (-cross(omega, params.J * omega) + M);
            
            dot_x(7:9, 1) = dot_angles;
            dot_x(10:12, 1) = dot_omega;
        case 'matrix'
            A = orthogonalize(x(7:9, 1:3)); % ортогонализация
            omega = x(10:12, 1);
            
            omegaMatrix = vec2matrix(omega);
            dotA = - omegaMatrix * A;
    
            dot_omega = inv(params.J) * (-cross(omega, params.J * omega) + M);

            dot_x(7:9, 1:3) = dotA;
            dot_x(10:12, 1) = dot_omega;
        case 'quat'
            quat = x(7:10, 1) / norm(x(7:10, 1)); % нормировка
            omega = x(11:13, 1);
            
            omega_quat = [0, omega];
            quat_dot = 1/2 * quatmultiply(quat, omega_quat);
            dot_omega = inv(params.J) * (-cross(omega, params.J * omega) + M); 
            
            dot_x(7:10, 1) = quat_dot;
            dot_x(11:13, 1) = dot_omega;
    end
end