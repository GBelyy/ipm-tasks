function [u, varargout] = calcAngularControl(r, v, quat, omega, t, params, sensors)
    usedJ = params.modelJ;

    switch lower(params.controlType)
        case 'bdot'
            mag = sensors.mag;
            B = mag.measure(t, [r;v;quat;omega], params);
            % bounding
            if t ~= 0
                dotB = (B - mag.prevValue) / params.dt;
            else
                dotB = [0;0;0];
            end

            mag.prevValue = B;
            u = -params.bdot_coef * cross(dotB, B);

            if (norm(params.bdot_coef * dotB) > params.torqueMoment)
                u = -cross(dotB, B) * params.torqueMoment / norm(dotB); %bounding by torque moment
            end

        case 'lyapunov'
            M_ext = calcExternalMoment([r; v; quat; omega], t, quat2dcm(quat'), params);
            
            [quatRef, omegaRef] = calcRefMotion(r, v, params);

            quatRel = quatmultiply(quatconj(quatRef'), quat')';
            q0 = quatRel(1);
            qVec = quatRel(2:4);

            A = quat2dcm(quatRel');
            omegaRel = omega - A * omegaRef;

            if t~= 0 
                dotOmegaRef = (omegaRef - params.omega_prev) / params.dt;
            else
                dotOmegaRef = params.omega_prev / params.dt;
            end

            u = - M_ext + cross(omega, usedJ * omega) + usedJ * A * dotOmegaRef ...
                - usedJ * cross(omegaRel, A * omegaRef) - 1/2 * params.kQ * q0 * params.modelJ * qVec...
                - params.kW * params.modelJ * omegaRel;
            
            varargout{1} = quatRel;
            varargout{2} = omegaRel;
            varargout{3} = omegaRef;
        case 'sliding'
            q0 = quat(1);
            qVec = quat(2:4);
            M_ext = calcExternalMoment([r; v; quat; omega], t, quat2dcm(quat'), params);
            u = -M_ext + cross(omega, usedJ * omega) - usedJ * params.P / params.lambda * (qVec + params.lambda * omega) - usedJ / (2 * params.lambda) *(q0 * omega + cross(qVec, omega));

        case 'magsliding'
            %https://www.keldysh.ru/papers/2014/prep2014_56_rus.pdf
            M_ext = calcExternalMoment([r; v; quat; omega], t, quat2dcm(quat'), params);

            mag = sensors.mag;
            B = mag.measure(t, [r;v;quat;omega], params);

            prevMat = params.lambda_prev;
            D = quat2dcm(quat');

            S = [D(2,3) - D(3,2);
                D(3,1) - D(1,3);
                D(1,2) - D(2,1)];
            dotS = [-D(1,3) * omega(3) + D(3,3) * omega(1) - D(1,2) * omega(2) + D(2,2) * omega(1);
                     D(1,1) * omega(2) - D(2,1) * omega(1) - D(2,3) * omega(3) + D(3,3) * omega(2);
                     D(2,2) * omega(3) - D(3,2) * omega(2) - D(3,1) * omega(1) + D(1,1) * omega(3)];

            a =(-params.dotLambda * usedJ * omega + params.lambda * (cross(omega, usedJ * omega) - M_ext) - ...
                prevMat * (usedJ * dotS + params.P * S) - params.lambda * params.P * omega) * params.dt + ...
                prevMat * usedJ * S;
            b = -usedJ * S;
            d = params.lambda * params.dt * B;

            e1 = d / norm(d);
            e3 = cross(d,b) / norm(cross(d,b));
            e2 = cross(e3,e1);

            T = [e1, e2, e3];

            a = T * a;
            b = T * b;
            d = T * d;

            currMat = zeros(3,3);
            currMat(1,1) = prevMat(1,1);
            currMat(1,2) = (-a(1) - currMat(1,1) * b(1))/b(2);
            currMat(2,1) = (-a(1) - currMat(1,1) * b(1))/b(2);
            currMat(2,2) = params.lambda0 + currMat(1,2)^2 / currMat(1,1);
            currMat(3,3) = prevMat(3,3);
            
            u = 1/(params.lambda * params.dt) * (-params.dotLambda * usedJ * omega + ...
                params.lambda * (cross(omega, usedJ * omega) - M_ext) - ...
                prevMat * (usedJ * dotS + params.P * S) - params.lambda * params.P * omega) * params.dt - ...
                currMat * usedJ * S + prevMat * usedJ * S;
        case 'no'
            u = [0;0;0];
    end
end