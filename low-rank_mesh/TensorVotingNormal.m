function [ tensorNormal ] = TensorVotingNormal( fNf, faceNormals, numF, init_iters, angle_thlocal, spatialW )
%TENSORVOTINGNORMAL Summary of this function goes here
%   Detailed explanation goes here
%faceNormals = gpuArray(faceNormals);

%num_neigh = size(fNf, 2);
Normal1 = faceNormals;

for iter=1:init_iters
    tensorNormal = zeros(numF, 3);
    parfor i=1:numF
        idxs = fNf{i,1};% should exclude the current point itself!!!
        %idxs = ballNeighs{i};% should exclude the current point itself!!!
        num_neigh = numel(idxs);
        tmpnormal = queryNormal(Normal1, idxs);
        Weight_local = exp( -(   (1-dot(repmat(Normal1(i,:),num_neigh,1), tmpnormal, 2)) ./ (1-cos(angle_thlocal))   ).^2 );
                                                                                                
                                                                                            %faceNormals(idxs(1:num_neigh),:)
        Weight_local = Weight_local .* spatialW{i,1}; %exp(-  sum( (repmat(Point.pos(i,:),num_neigh,1)-...
                                                           %          Point.pos(idxs(1:num_neigh),:)).^2, 2) ./ sigma_dis2 );
%     Tensor_local = zeros(3,3);
%     for j=1:num_neigh
%         Tensor_local = Tensor_local + Weight_local(j)* faceNormals(idxs(j),:)'*faceNormals(idxs(j),:);
%     end
    % same as the above
        %Tensor_local = faceNormals(idxs,:)' * spdiags(Weight_local,0,num_neigh,num_neigh) * faceNormals(idxs,:);
        %Tensor_local = tmpnormal' * spdiags(Weight_local,0,num_neigh,num_neigh) * tmpnormal;
        Tensor_local = tmpnormal' * diag(Weight_local,0) * tmpnormal;
        [eigenVec, eigenValues] = eig(Tensor_local);
        eigens = diag(eigenValues);
%         if(eigens<0)%debugging: to see if there are any negative values
%             fprintf('negative eigen values...%d\n', i);
%         end
        [~, eigenidx] = sort(eigens);% ascending order
        tensorNormal(i,:) = eigenVec(:,eigenidx(3))';
        %the eigenvector direction could be wrong: a flipped faceNormals is wrong
        consisDire = Normal1(i,:) * tensorNormal(i,:)';
        if consisDire < 0
            tensorNormal(i,:) = - tensorNormal(i,:);
        end
    end
    
    %tensorNormal = gather(tensorNormal);
    %if using this, essentially the combination of tensor filter and WNNM
    %uncomment when comparing with Bilateral Filtering
    
    %under the sense of multi iterations
    Normal1 = tensorNormal;
    
    
end


end

