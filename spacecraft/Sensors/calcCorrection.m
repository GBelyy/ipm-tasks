function [z, x_model, H] = calcCorrection(t, X, params, sensors)
    % observation matrises from https://www.keldysh.ru/papers/2014/prep2014_64.pdf
    switch (lower(params.filtertype))
        case 'ars+mag+sol'
            % measurements data
            ars = sensors.ars;
            mag = sensors.mag;
            sol = sensors.sol;
            D = quat2dcm(X(7:10)');

            B_meas = mag.measure(t, X, params);
            B_meas = B_meas / norm(B_meas);
            S_meas = sol.measure(t, X, params);
            S_meas = S_meas / norm(S_meas);
            W_meas = ars.measure(t, X, params);

            z = [B_meas; S_meas; W_meas];

            % model data
            B_i = calcMagneticField(X(1:6), t, params);
            B_model = D *  B_i/norm(B_i);
            JD = 2459580.5 + t/86400;
            S_i =  sun(JD)';
            S_model = D *  S_i/norm(S_i);
            W_model = X(11:13);

            x_model = [B_model; S_model; W_model];
            
            % obseravtion matrix
            H = [skewMat(B_model), zeros(3);
                 skewMat(S_model), zeros(3);
                 zeros(3), eye(3)];

        case 'star+ars'
            ars = sensors.ars;
            star = sensors.star;

            % measurements data
            quat_star = star.measure(t, X, params);
            Q_meas = quat_star(1:4);
            W_meas = ars.measure(t, X, params);
            
            z = [Q_meas; W_meas];
            
            % model data
            Q_model = X(7:10);
            W_model = X(11:13);
            x_model = [Q_model; W_model];
            
            % obseravtion matrix
            H = eye(6);
        otherwise
            error('Such type of sensor is not available.')
    end
end