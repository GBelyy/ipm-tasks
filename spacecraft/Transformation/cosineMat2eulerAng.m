function angles = cosineMat2eulerAng(cosineMat)
    angles = zeros(3, 1);

    e1 = cosineMat(1:3, 1);
    e2 = cosineMat(1:3, 2);
    e3 = cosineMat(1:3, 3);
    
    phi   = atan2(e1(3), e2(3));
    theta = acos(e3(3));
    psi   = atan2(e3(1), -e3(2));

    angles(1) = phi;
    angles(2) = theta;
    angles(3) = psi;
end