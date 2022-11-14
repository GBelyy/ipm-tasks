function kepElements = st2oe(stateVec, params)
    kepElements = zeros(6, 1);

    r = stateVec(1:3);
    v = stateVec(4:6);
    
    integrals = calcIntegrals(stateVec, params);
    h = integrals.h;
    c = integrals.c;
    f = integrals.f;

    n = c / norm(c); 
    l = cross([0; 0; 1], n) / norm(cross([0; 0; 1], n));
    j = cross(n, l);
    k = cross(n, f)/ norm(f);

    kepElements(1) = - params.earthGM/ (2 * h);         % semi major axis
    kepElements(2) = acos(n(3));                        % inclination
    kepElements(3) = atan2(l(2), l(1));                 % RAAN
    kepElements(4) = norm(f)/ params.earthGM;           % eccentricity
    kepElements(5) = atan2(dot(f, j), dot(f, l));       % arguement of perigee
    kepElements(6) = atan2(dot(r, k), dot(r, f));       % true anomaly

end

