function [M_ext, M_grav, M_atm, M_sol]  = calcExternalMoment(x, t, D, params)
    
    % geomagnetic field
    r_cck = D * x(1:3);
    M_grav = 3 * params.earthGM / norm(r_cck)^5 * cross(r_cck, params.J * r_cck);
    
    % resudial magnetic moment
    B_i = calcMagneticField(x(1:6), t, params);
    M_magn = cross(params.RMM, D * B_i);

    % athmosphere
    if isfield(params, 'atmOption') && params.atmOption
        [~, M_atm] = calcAtmResistance(x(1:3), x(4:6), D, t, params);
    else
        M_atm = [0;0;0];
    end
    % solar pressure
    if isfield(params, 'solarOption') && params.solarOption
        [~, M_sol] = calcSolarPressure(x(1:3), x(4:6), D, t, params);
    else
        M_sol = [0;0;0];
    end

    M_ext = M_grav + M_magn + M_atm + M_sol;
end