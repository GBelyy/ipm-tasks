function unitQuat = vec2Quat(vec)
    if vecnorm(vec) > 1
        unitQuat = [0; vec / norm(vec)];
    else
        unitQuat = [sqrt(1 - norm(vec)^2); vec];
    end
end