function [X_est, P_est]  = calcKalman(t1, t2, X, P_prev, params, sensors)
     % https://www.keldysh.ru/papers/2009/prep48/prep2009_48.pdf
    
     % state initializtion    
     t = linspace(t1, t2, 10);
     dt = abs(t(1) - t(2));
     N = length(t);

     X_pred = zeros(13, N);
     X_pred(:, 1) = X;

     % prediction 
     for i = 1:N-1 % prediction  of state to t2
         k1 = rhsAltitude(X_pred(:,i), t(i), params);
         k2 = rhsAltitude(X_pred(:,i) + dt / 2 * k1, t(i) + dt / 2, params);
         k3 = rhsAltitude(X_pred(:,i) + dt / 2 * k2, t(i) + dt / 2, params);
         k4 = rhsAltitude(X_pred(:,i) + dt * k3, t(i) + dt, params);
         X_pred(:,i + 1) = X_pred(:,i) + dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
         X_pred(7:10,i + 1) = X_pred(7:10,i + 1) / norm(X_pred(7:10,i + 1));         
     end
     rv = X_pred(1:6, N);

     x_pred = X_pred(7:13, N); % predicted angular state for t2 time step
     omega = x_pred(5:7);
     quat = x_pred(1:4);
     D = quat2dcm(quat');
     
     % evolution matrix from https://www.keldysh.ru/microsatellites/IAC_22_Ivanov.pdf  
     e = D * rv(1:3) / norm(rv(1:3));
     Fg = skewMat(e) * params.J * skewMat(e) - skewMat(params.J * e) * skewMat(e);
     Fx = skewMat(omega) * params.J - skewMat(params.J * omega);
     F = [    -skewMat(omega),      1/2 * eye(3);
          6 * params.earthGM / norm(rv(1:3))^5 * params.invJ * Fg,  (-params.invJ * Fx)];

     Phi = eye(6) + F * abs(t2 - t1);
     Q = params.Q;
     P_pred = Phi * P_prev * Phi' + Q;

     % correction
     [z, x_model, H] = calcCorrection(t2, X_pred(:, N), params, sensors);
     dw = z(5:7) - x_model(5:7);
     dQ = quatmultiply(quatconj(x_model(1:4)'), z(1:4)')';
     dq = dQ(2:4);

     dz = [dq; dw];
     % gain matrix
     K = P_pred * H'*inv(H * P_pred * H' + params.R); 

     x_corr = K * dz;
     q_corr = vec2Quat(x_corr(1:3));

     x_est = zeros(7, 1);
     x_est(1:4) = quatmultiply(quat', q_corr');
     x_est(5:7) = omega + x_corr(4:6);

     X_est = [rv; x_est];
     P_est = (eye(6) - K * H) * P_pred;
end