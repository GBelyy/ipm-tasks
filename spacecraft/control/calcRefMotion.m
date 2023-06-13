function [quatRef, omegaRef] = calcRefMotion(r, v, params)

    j3 = r/norm(r);
    j2 = cross(v, r) / norm(cross(v, r));
    j1 = cross(j2, j3) / norm(cross(j2,j3));

    coef = deltaFunc(isequal(params.gravOption, lower('j2')));
    J2coef = 3/2 * params.J2 * params.earthGM * params.earthRadius^2;
    a = coef * J2coef * r / norm(r)^5 * (5 * r(3)^2 / norm(r)^2 - 1) - 2 * coef / norm(r)^5 * [0; 0; r(3)];

    j3dot = (v - j3 * dot(v, j3)) / norm(r);
    hdot = cross(r, a);
    j2dot = (hdot - j2 * dot(hdot, j2)) / norm(cross(r,v));
    j1dot = cross(j2dot, j3) + cross(j2, j3dot);
    
    omegaRef = [dot(j2dot, j3);
                dot(j3dot, j1);
                dot(j1dot, j2)];
    Q = [j1'; j2'; j3'];

    quatRef = dcm2quat(Q)';
    if norm(quatRef - params.quat_prev) > norm(quatRef + params.quat_prev)
           quatRef = -quatRef;
    end
end