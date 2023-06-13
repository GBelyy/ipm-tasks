function [dot_x] = rhsAltitude(x, t, params, M_ctrl)
    if t == 8291
        disp(1)
    end
    if nargin < 4
        M_ctrl = [0;0;0];
    end
    D = quat2dcm(x(7:10, 1)');

    %if norm(x(1:3)) <= params.earthRadius
    %    disp(['Time of collision: ' num2str(t/3600) ' hours']);
    %    disp('Сollision with Earth is occured');
    %    quit
    %end

    if isfield(params, 'atmOption') && params.atmOption
        [F_atm, ~] = calcAtmResistance(x(1:3), x(4:6), D, t, params);
    else
        F_atm = [0;0;0];
    end
    % solar pressure
    if isfield(params, 'solarOption') && params.solarOption
        [F_sol, ~] = calcSolarPressure(x(1:3), x(4:6), D, t, params);
    else
        F_sol = [0;0;0];
    end
    F_perturb = F_atm + F_sol;
    a_perturb = F_perturb / params.mass; 

    [m, n] = size(x);
    dot_x = zeros(m, n);
    
    % Calculate orbital motion
    switch lower(params.gravOption)
        case '2bp'
            r = x(1:3, 1);
            v = x(4:6, 1);
            dot_x(1:3, 1) = v;
            dot_x(4:6, 1) = - params.earthGM / norm(r)^3 * r + a_perturb;
        case 'j2'
            r = x(1:3, 1);
            v = x(4:6, 1);
            coef = 3/2 * params.J2 * params.earthGM * params.earthRadius^2;
            accelerationJ2 = coef * r / norm(r)^5 * (5 * r(3)^2 / norm(r)^2 - 1) - 2 * coef / norm(r)^5 * [0; 0; r(3)];
            
            dot_x(1:3, 1) = v;
            dot_x(4:6, 1) = - params.earthGM / norm(r)^3 * r + accelerationJ2 + a_perturb;
        otherwise
            error('Wrong value gravOption');
    end
    
   switch lower(params.angParametrization)
       case 'quat'
            quat = x(7:10, 1); % нормировка
            omega = x(11:13, 1);
            
            M_ext = calcExternalMoment(x, t, D, params);
            
            M = M_ext + M_ctrl;
            omega_quat = [0; omega];
            quat_dot = (1/2 * quatmultiply(quat', omega_quat'))';
            dot_omega = params.invJ * (-cross(omega, params.J * omega) + M); 
            
            dot_x(7:10, 1) = quat_dot;
            dot_x(11:13, 1) = dot_omega;
        otherwise
            error('Wrong value angParametrization');
    end
end