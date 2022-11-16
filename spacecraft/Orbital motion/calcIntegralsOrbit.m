function integrals = calcIntegralsOrbit(stateVec, params)
    N = size(stateVec, 2);

    integrals = struct();
    h = zeros(1, N);
    c = zeros(1, N);
    f = zeros(1, N);
    potEnergy = zeros(1, N);

    for i = 1:N
        r = stateVec(1:3, i);
        v = stateVec(4:6, i);
        
        h(i) = dot(v, v) / 2 - params.earthGM / norm(r);
        c_vec = cross(r, v);
        c(i) = norm(c_vec);
        f(i) = norm(cross(v, c_vec) - params.earthGM * r/ norm(r));

        potEnergy(i) = params.earthGM / norm(r);
    end

    integrals.h = h;
    integrals.f = f;
    integrals.c = c;
    integrals.potEnergy = potEnergy;
end