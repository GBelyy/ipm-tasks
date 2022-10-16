function [dot_x] = rightSide(x, t, params)
    dot_x = zeros(6, 1);
    r = x(1:3);
    v = x(4:6);
    switch lower(params.option)
        case '2bp'
            dot_x(1:3) = v;
            dot_x(4:6) = - params.earthGM / norm(r)^3 * r;
        case 'j2'
            coef = 3/2 * params.J2 * params.earthGM * params.earthRadius^2;
            accelerationJ2 = coef * r / norm(r)^5 * (5 * r(3)^2 / norm(r)^2 - 1) - 2 * coef / norm(r)^5 * [0; 0; r(3)];
            
            dot_x(1:3) = v;
            dot_x(4:6) = - params.earthGM / norm(r)^3 * r + accelerationJ2;
    end
end