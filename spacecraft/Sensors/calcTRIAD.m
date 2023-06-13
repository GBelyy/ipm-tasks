function quat = calcTRIAD(t, X, params, mag, sol)
    %http://www.iki.rssi.ru/books/2013avanesov.pdf

    %ISK
    JD = 2459580.5 + t/86400;
    S_i =  sun(JD)';
    S = S_i/norm(S_i);
    rv = X(1:6);
    B_i = calcMagneticField(rv, t, params);
    B = B_i/norm(B_i);
    
    %SSK
    s = sol.measure(t, X, params);
    s = s/norm(s);
    b = mag.measure(t, X, params); 
    b = b/norm(b);

    C_I = cross(S,B);
    Gi = [S, B, C_I];     
    C_S = cross(s,b);
    Gs = [s,b,C_S];      

    D = Gs*Gi^(-1);
    
    crBS = cross(B,S);
    G = [B, crBS/norm(crBS), cross(B, crBS)/ norm(cross(B, crBS))];
    crbs = cross(b,s);
    g = [b, crbs/norm(crbs), cross(b, crbs)/ norm(cross(b, crbs))];
    
    A = g * G';
    A = orthogonalize(A);
    quat = dcm2quat(A)';
    
    if norm(quat - X(7:10)) > norm(quat + X(7:10))
           quat = -quat;
    end 
end