function [dot_x] = rightSide(x, t, params)
    switch lower(params.option)
        case 'simple'
            dot_x = [x(2); - params.g / params.l * sin(x(1))];
        case 'friction'
            dot_x = [x(2); -params.b * x(2) - params.g / params.l * sin(x(1))];
        case 'control'
            control = getControl(x, params);
            dot_x = [x(2); - params.g / params.l * sin(x(1)) + control];
    end
end