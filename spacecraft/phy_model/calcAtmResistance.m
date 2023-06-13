function [F, M]  = calcAtmResistance(r, v, D, t, params)
    normals = params.normals;

    eps = 0.1;
    alpha = 0.1;
    rho = calcAtmDensity(r, params);
    wEarth = [0;0;1] * params.omegaEarth; % angular velocity of Earth atmosphere

    F = zeros(3,1);
    M = zeros(3,1);
    
    vAtm = v - cross(wEarth, r);
    vAtm = vAtm * 1e3; % convert to m/c
    ev = vAtm/ norm(vAtm);
    
    for i = 1:6
        n_i = D' * normals(:,i);
        n_i = n_i/ norm(n_i);
        if dot(n_i, ev) > 0 % from cck to inertial
            R_i = 1/2 * params.overallSize * normals(:,i);
            S_i = params.overallSize * params.overallSize;
            F_i = - rho * dot(vAtm, vAtm) * dot(n_i, ev) *...
                  ((1-eps) * ev + ...
                    2 * eps * dot(ev, n_i) * n_i + ...
                   (1-eps) * alpha / norm(vAtm) * n_i) *  S_i;
            M_i = cross(R_i, D' * F_i);
        else
            F_i = zeros(3,1);
            M_i = zeros(3,1);
        end
        F = F + F_i;
        M = M + M_i;
    end
end