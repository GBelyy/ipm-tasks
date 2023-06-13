function changed_J = spoilTensor(J, sigma, koef)
    angle = [randn*sigma, randn*sigma, randn*sigma];
    quat = eul2quat(angle, 'XYZ');

    rot_mat = quat2dcm(quat);
    coef_Mat = randn(3,1) * koef + eye(3);
    changed_J = J .* coef_Mat;

    changed_J = rot_mat' * changed_J * rot_mat;
end