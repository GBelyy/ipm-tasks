function integrals = calcIntegralsOrbit(stateVec, params)
    N = size(stateVec, 2);

    integrals = struct();

    r = stateVec(1:3);
    v = stateVec(4:6);

    h = dot(v, v) / 2 - params.earthGM / norm(r);
    c_vec = cross(r, v);
    c = norm(c_vec);
    f_vec = cross(v, c_vec) - params.earthGM * r/ norm(r);
    f = norm(f_vec);

    potEnergy = params.earthGM / norm(r);

    integrals.h = h;
    integrals.f = f;
    integrals.c = c;
    integrals.c_vec = c_vec;
    integrals.f_vec = f_vec;
    integrals.potEnergy = potEnergy;
end