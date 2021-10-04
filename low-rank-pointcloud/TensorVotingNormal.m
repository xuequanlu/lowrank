function [ tensorNormal ] = TensorVotingNormal( localNNpoints, Normal, numPoints, init_iters, angle_thlocal, spatialW )
%TENSORVOTINGNORMAL Summary of this function goes here
%   Detailed explanation goes here

%num_neigh = size(localNNpoints, 2);
Normal1 = Normal;

for iter=1:init_iters
    tensorNormal = zeros(numPoints, 3);
    parfor i=1:numPoints
        idxs = localNNpoints(i, :);
        
        num_neigh = numel(idxs);
        tmpnormal = queryNormal(Normal1, idxs);
        Weight_local = exp( -(   (1-dot(repmat(Normal1(i,:),num_neigh,1), tmpnormal, 2)) ./ (1-cos(angle_thlocal))   ).^2 );
                                                                                                
        Weight_local = Weight_local .* spatialW(:,i);

        Tensor_local = tmpnormal' * diag(Weight_local,0) * tmpnormal;
        [eigenVec, eigenValues] = eig(Tensor_local);
        eigens = diag(eigenValues);

        [~, eigenidx] = sort(eigens);% ascending order
        tensorNormal(i,:) = eigenVec(:,eigenidx(3))';
        %the eigenvector direction could be wrong: a flipped normal is wrong
        consisDire = Normal1(i,:) * tensorNormal(i,:)';
        if consisDire < 0
            tensorNormal(i,:) = - tensorNormal(i,:);
        end
    end
    
    Normal1 = tensorNormal;
    
    
end


end

