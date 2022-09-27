function [control] = getControl(x, params)
% Function returns control action depending on system option
    if strcmp(params.option, 'control')
        control = params.g / params.l * sin(x(1)) ...
                  - params.koefPhi * (x(1) - params.xRef(1)) ...
                  - params.koefOmega * (x(2) - params.xRef(2));
    else
        control = 0;
    end
end