function [dot_x] = rhsAngle(x, t, params)
    [m, n] = size(x);
    dot_x = zeros(m, n);
    switch lower(params.angParametriztion)
        case 'euler'
            phi = x(1, 1);
            theta = x(2, 1);
            psi = x(3, 1);
            omega = x(4:6);

            dot_angles(1) = omega(1) * cos(phi) - omega(2) * sin(phi);
            dot_angles(2) = (omega(1) * sin(phi) + omega(2) * cos(phi)) / sin(theta);
            dot_angles(3) = omega(3) - cot(theta) * (omega(1) * sin(phi) + omega(2) * cos(phi));

            dot_omega = inv(params.J) * (-cross(omega, params.J * omega) + M);
            
            dot_x(1:3, 1) = dot_angles;
            dot_x(4:6, 1) = dot_omega;
        case 'matrix'
            A = orthogonalize(x(1:3, 1:3)); % ортогонализация
            omega = x(4:6, 1);
            
            omegaMatrix = vec2matrix(omega);
            dotA = - omegaMatrix * A;
    
            dot_omega = inv(params.J) * (-cross(omega, params.J * omega) + M);

            dot_x(1:3, 1:3) = dotA;
            dot_x(4:6, 1) = dot_omega;
        case 'quat'
            quat = x(1:4, 1) / norm(x(1:4, 1)); % нормировка
            omega = x(5:7, 1);
            
            omega_quat = [0, omega];
            quat_dot = 1/2 * quatmultiply(quat, omega_quat);
            dot_omega = inv(params.J) * (-cross(omega, params.J * omega) + M); 
            
            dot_x(1:4, 1) = quat_dot;
            dot_x(5:7, 1) = dot_omega;
    end
end