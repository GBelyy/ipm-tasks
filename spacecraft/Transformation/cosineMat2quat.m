function quat = cosineMat2quat(cosineMat)
    quat = zeros(4, 1);

    e1 = cosineMat(1:3, 1);
    e2 = cosineMat(1:3, 2);
    e3 = cosineMat(1:3, 3);

    quat(1) = 1/2 * sqrt(1 + e1(1) + e2(2) + e3(3));
    quat(2) = (e2(3) - e3(2)) / (4 * quat(1));
    quat(3) = (e3(1) - e1(3)) / (4 * quat(1));
    quat(4) = (e1(2) - e2(1)) / (4 * quat(1));

    % случай малой скалярной части
    if quat(1) <= 1e-6
        if quat(2) >= 1e-6
            quat(2) = 1/2 * sqrt(1 + e1(1) - e2(2) - e3(3));
            quat(1) = (e2(3) - e3(2)) / (4 * quat(2));
            quat(3) = (e1(2) + e2(1)) / (4 * quat(2));
            quat(4) = (e1(3) + e3(1)) / (4 * quat(2));

        elseif quat(3) >= 1e-6
            quat(3) = 1/2 * sqrt(1 - e1(1) + e2(2) - e3(3));
            quat(1) = (e3(1) - e1(3)) / (4 * quat(3));
            quat(2) = (e1(2) + e2(1)) / (4 * quat(3));
            quat(4) = (e2(3) + e3(2)) / (4 * quat(3));
        elseif quat(4) >= 1e-6
            quat(4) = 1/2 * sqrt(1 - e1(1) - e2(2) + e3(3));
            quat(1) = (e1(2) - e2(1)) / (4 * quat(4));
            quat(2) = (e1(3) + e3(1)) / (4 * quat(4));
            quat(3) = (e2(3) + e3(2)) / (4 * quat(4));
        end
    end
end