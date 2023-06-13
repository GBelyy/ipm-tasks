function xRTN = eci2rtn(x)
    r = x(1:3);
    v = x(4:6);

    e3 = r/norm(r);
    e2 = cross(r,v) / norm(cross(r,v));
    e1 = cross(e2, e3);

    D = [e1, e2, e3]; % ECI to RTN matrix;
    xRTN = D * x;
end