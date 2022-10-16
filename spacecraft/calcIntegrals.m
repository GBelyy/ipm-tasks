function integrals = calcIntegrals(stateVec, params)

    r = stateVec(1:3);
    v = stateVec(4:6);
    
    integrals = struct();
    
    integrals.h = dot(v, v) / 2 - params.earthGM / norm(r);
    integrals.c = cross(r, v);
    integrals.f = cross(v, integrals.c) - params.earthGM * r/ norm(r);
end