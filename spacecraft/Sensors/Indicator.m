classdef Indicator < handle
    properties
        indicatorType
        prevValue
        params
    end

    methods
        function this = Indicator(type,params)
            this.indicatorType = lower(type);
            this.params = params;
        end

        function estimate = measure(this, t, X, params)
            orientReal = this.params.orientReal;
            orientModel = this.params.orientModel;
            E = this.params.Expect;
            sigma = this.params.sigma;
            
            error = zeros(3,1);
            for i = 1:3
                vec = randn;
                error(i) = randn*sigma + E; 
            end
            quat = X(7:10);
            D = quat2dcm(quat');
            switch(this.indicatorType)
                case 'ars'
                    omega = X(11:13);
                    omega_err = omega + error;
                    estimate = orientModel' * orientReal * omega_err;
                case 'magnetometer'
                    rv = X(1:6);
                    B_i = calcMagneticField(rv, t, params);
                    B_cck = D * B_i;

                    estimate = orientModel' * (orientReal' * B_cck + error);
                case 'sunsensor'
                    JD = 2459580.5 + t/86400;

                    e3 = sun(JD)/norm(sun(JD));
                    if norm(cross([1;0;0],e3)) < 1e-12
                        e1 = cross([1;0;0],e3);
                    else
                        e1 = cross([0;1;0],e3);
                    end
                    
                    e1 = e1/(norm(e1));
                    e2 = cross(e3,e1);
                    C = [e1;e2;e3];
                    theta = randn * sigma + E;
                    phi = pi*rand;
                    sol_err = [sin(theta)*cos(phi);
                               sin(theta)*sin(phi);
                               cos(theta)];

                    estimate = orientModel'* orientReal * D * C' * sol_err;
                case 'star'
                    % TBD - right quat noising
                    star_err = quatmultiply(dcm2quat(orientReal), vec2Quat(error)')';
                    spoiledQuat = quatmultiply(quat', star_err')';
                    estimate = quatmultiply(spoiledQuat', quatconj(dcm2quat(orientModel)))';
                otherwise
                    warning('Sensor type not founded.')
                    estimate = [];
            end
        end
    end
end