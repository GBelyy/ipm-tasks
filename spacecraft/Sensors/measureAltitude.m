function [quat, omega] = measureAltitude(t, X, params, sensors, varargin)
    switch lower(params.measModel)
        case 'ideal'
            quat = X(7:10);
            omega = X(11:13);
        case 'triad'
            omega = X(11:13);

            mag = sensors.mag;
            sol = sensors.sol;

            quat = calcTRIAD(t, X, params, mag, sol);
        case 'kalman'
            % TODO: realize getVararginParam

            % sequence of input varargin:
            % varargin{1} - t1 - time moment of previous calculation
            % varargin{2} - P_prev - covariation matrix from previous time step
            
            t1 = varargin{1};
            P_prev = varargin{2};
            warning('Method in developing ...');
            
            [X_est, ~]  = calcKalman(t1, t2, X, P_prev, params, sensors);
            quat = X_est(7:10);
            omega = X_est(11:13);
    end
end