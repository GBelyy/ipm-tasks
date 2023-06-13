function [F, M]  = calcSolarPressure(r, v, D, t, params)
    % https://www.keldysh.ru/papers/2017/prep2017_78_rus.pdf
    
    % solar panels are located on the surfaces of spacecraft and has area of these surfaces 
    normals = params.panelNormals;
    panelArea = params.panelArea;
    
    c = 299793.1 * 1e3; %speed of light, km/c
    F_s = 1367; % solar constant, W/m2

    alpha = 0.1;
    mu = 0.5;
    
    JD = 2459580.5 + t/86400;
    sunDir = sun(JD);
    sunDir = sunDir'/norm(sunDir);

    r_s = sunDir; % sun directed vector of spacecraft's center of mass
    r_s = r_s / norm(r_s);

    F = zeros(3,1);
    M = zeros(3,1);
    for i = 1:size(normals, 2)
        n_i = D' * normals(:,i);
        n_i = n_i/ norm(n_i);
        if dot(r_s, n_i) > cos(params.thetaMax)
            R_i = 1/2 * params.overallSize * normals(:,i);
            S_i = panelArea;
            F_i = - F_s/c * dot(r_s, n_i) * ...
                  ((1-alpha) * r_s + ...
                  2 * alpha * mu * dot(r_s, n_i) * n_i + ...
                  alpha * (1 - mu) * (r_s + 2/3 * n_i)) * ...
                  S_i;
            M_i = cross(R_i,D' * F_i);
        else
            F_i = zeros(3,1);
            M_i = zeros(3,1);
        end
        F = F + F_i;
        M = M + M_i;
    end
end