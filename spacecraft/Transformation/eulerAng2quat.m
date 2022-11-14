function quat = eulerAng2quat(angles)
    quat = zeros(4, 1);
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);

    quat(1) = cos(theta/2) * cos((psi + phi) / 2);
    quat(2) = sin(theta/2) * cos((psi - phi) / 2);
    quat(3) = sin(theta/2) * sin((psi - phi) / 2);
    quat(4) = cos(theta/2) * sin((psi + phi) / 2);
end