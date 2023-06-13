function B_inertial = calcMagneticField(rv, t, params)
    % https://www.keldysh.ru/microsatellites/Bachelor_Thesis_Pichuzhkina.pdf
    r = rv(1:3);

    lambda2 = params.omegaEarth * t;
    delta1 = deg2rad(params.theta);
    unitVec = r/norm(r);
    vec = unitVec(1) * sin(delta1) * sin(lambda2) - unitVec(2) * sin(delta1) * cos(lambda2) + unitVec(3) * cos(delta1);
    
    B_inertial = params.magMuEarth / norm(r)^3 * [sin(lambda2) * sin(delta1) - 3 * vec * unitVec(1);
                                                  -cos(lambda2) * sin(delta1) - 3 * vec * unitVec(2);
                                                  cos(delta1) - 3 * vec * unitVec(3)];
end