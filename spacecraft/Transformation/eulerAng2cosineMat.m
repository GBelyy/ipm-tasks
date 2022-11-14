function cosineMat = eulerAng2cosineMat(angles)
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);

    e1 = [cos(psi) * cos(phi) - sin(psi) * cos(theta) * sin(phi);
          sin(psi) * cos(phi) + cos(psi) * cos(theta) * sin(phi);
          sin(theta) * sin(phi)];
    e2 = [-cos(psi) * sin(phi) - sin(psi) * cos(theta) * cos(phi);
          -sin(psi) * sin(phi) + cos(psi) * cos(theta) * cos(phi);
          sin(theta) * cos(phi)];
    e3 = [sin(psi) * sin(theta);
          -cos(psi) * sin(theta);
          cos(theta)];

    cosineMat = [e1, e2, e3];
end