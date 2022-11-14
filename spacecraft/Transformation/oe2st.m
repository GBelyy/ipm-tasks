function stateVec = oe2st(kepElements, params)
    stateVec = zeros(6, 1);

    a = kepElements(1);
    inc = kepElements(2);
    RAAN = kepElements(3);
    e = kepElements(4);
    omega = kepElements(5);
    trueAnomaly = kepElements(6);
    
    p = a * (1 - e ^ 2);
    r = p / (1 + e * cos(trueAnomaly));
    
    % components of state vector in perifocal coordinate system
    r_pqw = [r * cos(trueAnomaly); ...
             r * sin(trueAnomaly); ...
             0];

    v_pqw = [-sqrt(params.earthGM / p) * sin(trueAnomaly); ...
             sqrt(params.earthGM / p) * (e + cos(trueAnomaly)); ...
             0];

    rotationMtx = rotMatrix(RAAN, 3) * rotMatrix(inc, 1) * rotMatrix(omega, 3);

    r = rotationMtx * r_pqw;
    v = rotationMtx * v_pqw;
    
    stateVec(1:3) = r;
    stateVec(4:6) = v;
end

function matrix = rotMatrix(alpha, axis)
    switch axis
        case 1
            matrix = [1, 0, 0; ...
                      0, cos(alpha), -sin(alpha); ...
                      0, sin(alpha), cos(alpha)];
        case 2
            matrix = [cos(alpha), 0, sin(alpha); ...
                      0, 1 0; ...
                      -sin(alpha), 0, cos(alpha)];
        case 3
            matrix = [cos(alpha), -sin(alpha), 0; ...
                      sin(alpha), cos(alpha), 0; ...
                      0, 0, 1];
        otherwise
    end
end