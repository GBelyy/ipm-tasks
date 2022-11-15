function ortMat = orthogonalize(mat)
    [m, n] = size(mat);
    ortMat=zeros(m, n); 
    
    ortMat(:, 1)=mat(:, 1);
    for i=2:n
        for j=1:i-1
            ortMat(:, i) = ortMat(:, i) - dot(mat(:, i),ortMat(:, j)) / dot(ortMat(:, j),ortMat(:, j)) * ortMat(:, j);
        end
        ortMat(:, i) = ortMat(:, i) + mat(:, i);
    end
    
    % normalization
    for k=1:n
         ortMat(:,k)=ortMat(:,k)/norm(ortMat(:,k));
    end
end