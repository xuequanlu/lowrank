function [ spatialW ] = computeSpatialW(localNNpoints, Position, numPoints, sigma_dis2 )
   
    num_neigh = size(localNNpoints, 2);%the same for each point
    spatialW = zeros(num_neigh, numPoints);
    
    parfor i=1:numPoints
        idxs = localNNpoints(i, :);
        num_neigh = numel(idxs);
        tmpPos = queryPos(Position, idxs);
        weight = exp(-  sum( (repmat(Position(i,:),num_neigh,1) - tmpPos).^2, 2) ./ sigma_dis2(i) );
        spatialW(:, i) = weight;
        
    end
    
    
end

