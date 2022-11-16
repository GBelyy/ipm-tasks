function integrals = calcIntegralsAngle(stateVec, params)
    sizes = size(stateVec);
    N = sizes(end);
    
    integrals = struct();
    K = zeros(1, N);
    energy = zeros(1, N);
    k = int8(isequal(params.angParametrization,'quat')) + 4;
    switch lower(params.integralsEvent)
        case 'euler'
            for i = 1:N
                omega = stateVec(k:k+2, 1, i);
                K_vec = params.J * omega;
                K(i) = norm(K_vec);
                energy(i) = dot(omega, K_vec) / 2;
            end
            
            integrals.K = K;
            integrals.energy = energy;
        case 'lagrange'
            for i = 1:N
                omega = stateVec(k:k+2, 1, i);
                angles = transform2angles(stateVec(1:k - 1, :, i), params);
                rotateMat = inv(eulerAng2cosineMat(angles));

                K_vecFCS = params.J * omega;     % fixed coordinate system
                K_vecICS = rotateMat * K_vecFCS; % inertial coordinate system
                Ke3(i) = K_vecFCS(3);
                Kz(i) = K_vecICS(3);
                energy(i) = dot(omega, K_vecFCS) / 2 + params.m * params.g * params.l * cos(angles(2));
            end
            
            integrals.Ke3 = Ke3;
            integrals.Kz = Kz;
            integrals.energy = energy;
    end
    end

function angles = transform2angles(state, params)
    switch lower(params.angParametrization)
        case 'euler'
            angles = state;
        case 'cosinemat'
            angles = cosineMat2eulerAng(state);
        case 'quat'
            angles = quat2eulerAng(state);
        otherwise
            error('Неправильное имя параметра angParametrization');
    end
end