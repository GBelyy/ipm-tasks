function angleRad = dms2rad(dms)
    grads = dms(1);
    minutes = dms(2);
    secs = dms(2);
    
    angleRad = (grads + minutes/60 + secs/3600) * pi/180;
end