function [centroids] = computeCentroidofTriangles(V,F)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    [numF, dim] = size(F);
    centroids = zeros(numF, 3);
    for i=1:numF
        c = V(F(i,1),:) + V(F(i,2),:) + V(F(i,3),:);
        c = c / 3.0;
        centroids(i,:) = c;
    end

end

