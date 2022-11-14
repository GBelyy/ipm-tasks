function angles = quat2eulerAng(q)
    angles = zeros(3, 1);
    
    lamPos = atan2(q(4), q(1));
    lamNeg = atan2(q(3), q(2));

    phi   = lamPos - lamNeg;
    theta = 2 * atan(sqrt((q(2)^2 + q(3)^2) / (q(1)^2 + q(4)^2)));
    psi   = lamPos + lamNeg;

    angles(1) = phi;
    angles(2) = theta;
    angles(3) = psi;
end