function [dot_x] = rhsAngle(x, t, params)
    [m, n] = size(x);
    dot_x = zeros(m, n);
    
    switch lower(params.angParametrization)
        case 'euler'
            M = (cross([0;0;1], eulerAng2cosineMat(x(1:3,1)) * [0;0;-1]) * params.m * params.g * params.l) * double(isequal(params.integralsEvent,'lagrange'));
            % angular velocity in
            phi = x(1, 1);
            theta = max(x(2, 1), 1e-5);
            psi = x(3, 1);
    
            omega = x(4:6, 1);
            dot_angles = zeros(3,1);
    
            dot_angles(1) = omega(3) - cot(theta) * (omega(1) * sin(phi) + omega(2) * cos(phi));
            dot_angles(2) = omega(1) * cos(phi) - omega(2) * sin(phi);
            dot_angles(3) = (omega(1) * sin(phi) + omega(2) * cos(phi)) / sin(theta);
    
            dot_omega = params.invJ * (-cross(omega, params.J * omega) + M);
    
            dot_x(1:3, 1) = dot_angles;
            dot_x(4:6, 1) = dot_omega;
        case 'cosinemat'
            A = orthogonalize([x(1:3), x(4:6), x(7:9)]); % ortogonolization
            M = (cross([0;0;1],A * [0;0;-1]) * params.m * params.g * params.l) * double(isequal(params.integralsEvent,'lagrange'));
            omega = x(10:12, 1);
    
            omegaMatrix = vec2matrix(omega);
            dotA = - omegaMatrix * A;
    
            dot_omega = params.invJ * (-cross(omega, params.J * omega) + M);
    
            dot_x(1:9) = [dotA(:, 1); dotA(:, 2); dotA(:, 3)];
            dot_x(10:12, 1) = dot_omega;
        case 'quat'
            quat = x(1:4, 1)/ norm(x(1:4, 1)); % normalization
            M = (cross([0;0;1], quat2dcm(quat') * [0;0;-1]) * params.m * params.g * params.l) * double(isequal(params.integralsEvent,'lagrange'));
            omega = x(5:7, 1);
            omega_quat = [0; omega];
            quat_dot = (1/2 * quatmultiply(quat', omega_quat'))';
            dot_omega = params.invJ * (-cross(omega, params.J * omega) + M);
    
            dot_x(1:4, 1) = quat_dot;
            dot_x(5:7, 1) = dot_omega;
        otherwise
            error('Неправильное имя параметра angParametrization');
    end
end