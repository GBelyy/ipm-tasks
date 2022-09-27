function [dot_x] = rightSide(x, t, params, control)
% Function returns right side of motion equataion of mathematical pendulum
% depending on system option
    switch lower(params.option)
        case 'simple'
            dot_x = [x(2); - params.g / params.l * sin(x(1))];
        case 'friction'
            dot_x = [x(2); -params.b * x(2) - params.g / params.l * sin(x(1))];
        case 'control'
            dot_x = [x(2); - params.g / params.l * sin(x(1)) + control];
        otherwise
            error([params.option 'is not accepted. Choose one from simple/friction/control.'])
    end
end