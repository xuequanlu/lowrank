% for mesh denoising
%face faceNormals smoothing
%%
[V, F] = readPLY3('data/dod_mesh+gmm(0.3)+pre.ply');
F = F +1;%index start from 1

[numV, colV] = size(V);
[numF, colF] = size(F);

%color = ones(numF,4)*0;

%use embedded functions instead
V_col1 = V(F(:,1), :);
V_col2 = V(F(:,2), :);
V_col3 = V(F(:,3), :);

centroidPos = (V_col1 + V_col2 + V_col3) ./ 3.0;
faceNormals = cross( (V_col2-V_col1), (V_col3-V_col2), 2 );
faceNorms = repmat( sqrt( sum(faceNormals.^2, 2) ), 1, 3);
faceNormals = faceNormals ./ faceNorms;

tree = KDTreeSearcher( centroidPos );% build the tree of faces (centroid)

% 
angle_thnonlocal = 30.0 %
angle_minthnonlocal = 15.0
reduce_thnonlocal = 1.0;
cos_thnonlocal = cos(30.0*pi/180.0);
cos_maxthnonlocal = cos(15.0*pi/180.0);

angle_thlocal1 = 30.0 
angle_minthlocal1 = 15.0
reduce_thlocal1 = 1.0;
cos_thlocal1 = cos(30.0*pi/180.0);
cos_maxthlocal1 = cos(15.0*pi/180.0);

reduce_rate = 1.1;

angle_thlocal = 30.0 /180.0*pi;%for Tensor Voting: fixed
Patchsize_th = 1e10;
%a large enough region, can be tuned by users: sugget large number like 150 for CAD-like models, small value like 20 for models with details
knn = 150 

% Knn
[neighborPoints, dis_knn] = knnsearch(tree, centroidPos, 'K', knn);%fast: increasing distance order
%1-ring local structure
% vertex-to-face neighbors
vNf = cell(numV,1);
AdjFV = sparse(numF, numV);
for i=1:numF
    c1 = F(i,1);
    c2 = F(i,2);
    c3 = F(i,3);
    AdjFV(i, c1) = 1.0;
	AdjFV(i, c2) = 1.0;
	AdjFV(i, c3) = 1.0;
end
for i=1:numV
    [row,col] = find( AdjFV(:,i) );
    vNf{i,1} = row; % column
end
clear AdjFV;

% you can choose to use 1-ring or 2-ring as local structure
% 1-ring face-to-face neighbors
fNf1 = cell(numF, 1);
for i=1:numF
    c1 = F(i,1);
    c2 = F(i,2);
    c3 = F(i,3);
    %1-ring face neighbors
    fNf_i = [vNf{c1,1}; vNf{c2,1}; vNf{c3,1}];
    fNf_i = unique(fNf_i, 'stable');% should be continuous
    fNf_i(fNf_i==i) = [ ];% remove the current face
    fNf1{i,1} = fNf_i; % column vector
    
    %1-ring neighbors as local structure
    %fNf1{i,1} = fNf_i'; % row vector
end
%fNf = fNf1;
% %2-ring face neighbors, roughly about 30
fNf2 = cell(numF, 1);
for i=1:numF
    fNf_i = fNf1{i,1};
    fNf2_i = [fNf_i; cell2mat( fNf1(fNf_i) )] ;%2-ring
    fNf2_i = unique( fNf2_i, 'stable' );
    fNf2_i(fNf2_i==i) = [ ]; %remove current face
    fNf2{i,1} = fNf2_i' ;% comlun to row vector
end
fNf = fNf2;
clear fNf2;



sigma_dis2 = zeros(numF,1);
for i=1:numF
    num_fNf = numel( fNf{i,1} );
    dis2_1ring = sum (repmat(centroidPos(i,:), num_fNf, 1) - centroidPos(fNf{i,1}, :).^2, 2 );
    dis2_1ring = sort(dis2_1ring);%ascending order
    sigma_dis2(i,1) = 4.0* dis2_1ring(num_fNf);% the largest distance
end

%% (1) initialization for normals: with anisotropic Tensor voting (optional)
fprintf('computing spatialW...');
tic;
spatialW = computeSpatialW(fNf1, centroidPos, numF, sigma_dis2); % only compute once, using 1-ring could be enough
toc;
init_iters = 1;%typically 1 iteration, since it is initialization
fprintf('computing TensorVotingNormal for initial...');
tic;%using 1-ring could be enough
faceNormals = TensorVotingNormal( fNf1, faceNormals, numF, init_iters, angle_thlocal, spatialW );
toc;

allFactors = importdata('decomposedInteger.mat');

localIso = cell(1, numF);
simPatches = cell(1, numF);
%%
max_outerIter = 5;
outerIter = 0;

tic;
while (outerIter<max_outerIter)
%% (2) find similar patches for each point
%local representative faceNormals of each patch, not the faceNormals of this point!
   init_iters = 1;
   localstru = TensorVotingNormal( fNf1, faceNormals, numF, init_iters, angle_thlocal, spatialW );

parfor i = 1:numF
    idxs = neighborPoints(i, :); %enough large 'window'
    numofNeigh = numel(idxs);
    tmp_localstru = queryLocalstru(localstru, idxs);
    
    % record the isotropic neighbors consistent with localstructure
    idxs_ball = fNf{i,1};
    numofballNeigh = numel(idxs_ball);
    tmp_normal = queryNormal(faceNormals, idxs_ball);
    
    inAngle =  dot( repmat(localstru(i,:),numofballNeigh,1), tmp_normal, 2);
    
    [inAngle,idss] = sort(inAngle,'descend'); %decreasing order for cos(), increasing order in angle
    inAng1 = inAngle(inAngle >= cos_thlocal1);
    localIso{i} = [i, idxs_ball(:, idss(1:numel(inAng1),:) ) ];%should include the current point itself

    % find similar non-local patches (via a large enough 'window')
    
    interAngle = dot( repmat(localstru(i,:),numofNeigh,1), tmp_localstru, 2) ;
    [interAngle, idxx] = sort(interAngle, 'descend');%increasing order in angle
    interAngle1 = interAngle(interAngle >= cos_thnonlocal);%angle_thnonlocal
    simPatches{i} = idxs( :, idxx(1:numel(interAngle1),:) ); %increasing error order, include itself
end

%% (3) for each point, perform nuclear norm minimization (low rank matrix approximation)

sumNormals = zeros(numF,3);
sumWeights = zeros(numF,1);

%record the abnormal point ids
abpoints = sparse(numF,1);
ccc = 0;
parfor i=1:numF
    idxs = simPatches{i};
    nonlocal_patchsize = numel(idxs);
    
    % acquire the similar patches for the i-th point
    tmpLocalIso = queryLocalIso( localIso, idxs(1,1:nonlocal_patchsize) );
    vec_tmpoints = cell2mat( tmpLocalIso );
    
    num_actual = numel(vec_tmpoints);
    
    
    if num_actual <=5 %for debugging
        fprintf('Iter: %d, too little...%d, %d\n', outerIter, i, num_actual);
    end
    if num_actual <=5 && outerIter>=2
        %handle abnormal points
        num_fNf1 = numel( fNf1{i,1} );
        curr_1ring = fNf1{i,1};
        for tt=1:num_fNf1% search 1-ring neighbors
            tmInd = curr_1ring(tt);
            if ~any(tmInd==vec_tmpoints)
                abpoints(i,1) = tmInd;
                break;
            end
        end
    end
    
    [rowNum_patch, colNum_patch, Patchsize] = decomposeInteger( num_actual*3, allFactors );
    vec_tmpoints = vec_tmpoints(:,1:Patchsize);

    matVn_patch = queryNormal(faceNormals, vec_tmpoints);%faceNormals(vec_tmpoints, :);
    num_actual = numel(vec_tmpoints);
    
    
    matVn_patchX = matVn_patch;
    
    %xxxx...yyyy...zzzz...
    matVn_patchX = reshape(matVn_patchX, rowNum_patch, colNum_patch);
   
    
    % denoise each patch
    iter = 1;% for mesh and scanned models, 1 iter is enough
    whichIter = 0;
    while whichIter < iter
       
        [U,SigmaX,V] = svd(matVn_patchX, 'econ');
        
        for kk=1:3
            if kk==1
                Temp = diag(SigmaX);
            end
            W_VecX = exp(- (2.0*Temp./Temp(1)).^2 ); 
            SigmaX1 = soft(SigmaX, diag(W_VecX));
            Temp = diag(SigmaX1);
        end
    
        X =  U*SigmaX1*V';
        
        % standardization for normals
        
        patch_n = reshape(X, [], 3);% N * 3
        patchnorms = repmat( sqrt(sum(patch_n.^2, 2)), 1, 3 );
        patch_n = patch_n ./ patchnorms;
        %
        matVn_patchX = reshape(patch_n, rowNum_patch, colNum_patch);%N * M
        
        whichIter = whichIter+1;
    end

    
    %
    vec_true = unique(vec_tmpoints);
    realsize = numel(vec_true);
    reco_id = zeros(realsize, 1);
    reco_n = zeros(realsize, 3);
    for tt=1:realsize
        
        tt2 = (vec_true(tt)==vec_tmpoints);%logic index
        reco_id(tt,1) = sum(tt2);
        
        reco_n(tt,:) = sum( patch_n(tt2, :), 1);
    end
    %
    tmpPatchWeights = zeros(numF, 1);
    tmpPatchWeights(vec_true,:) = tmpPatchWeights(vec_true,:) + reco_id;
    sumWeights = sumWeights + tmpPatchWeights;
    %
    tmpPatchNormal = zeros(numF, 3);
    tmpPatchNormal(vec_true,:) = tmpPatchNormal(vec_true,:) + reco_n;
    sumNormals = sumNormals + tmpPatchNormal;
    
end



normal_preIter = faceNormals;%record the original normals
% average and standardization at each outer iteration
faceNormals = sumNormals ./ ( repmat(sumWeights,1,3) );% average, no influences
Allnorms = repmat( sqrt(sum(faceNormals.^2, 2)), 1, 3 ) ;
faceNormals = faceNormals ./ Allnorms;


%deal with specific abnormal points, after computing new normals
nnz_abpoints = nnz(abpoints);
if nnz_abpoints >= 1
    [rowc, colc, valuec] = find(abpoints);
    tmp = faceNormals(valuec, :);
    faceNormals(rowc, :) = tmp;
end

%% reduce angle thresholds iteratively
fprintf('Iter: %d, angle_thnonlocal: %d; angle_thlocal1: %d\n', outerIter, angle_thnonlocal, angle_thlocal1);
fprintf('Iter: %d, cos_thnonlocal: %d; cos_thlocal1: %d\n', outerIter, cos_thnonlocal, cos_thlocal1);


%
reduce_thnonlocal = reduce_thnonlocal * reduce_rate;
angle_thnonlocal = angle_thnonlocal - reduce_thnonlocal;
if angle_thnonlocal > angle_minthnonlocal
    cos_thnonlocal = cos( (angle_thnonlocal)*pi/180 );
else
    angle_thnonlocal = angle_minthnonlocal;
    cos_thnonlocal = cos_maxthnonlocal;
end
%use cos() instead of angles, to reduce computation
 reduce_thlocal1 = reduce_thlocal1 * reduce_rate;
 angle_thlocal1 = angle_thlocal1 - reduce_thlocal1;
if angle_thlocal1 > angle_minthlocal1
    cos_thlocal1 = cos( (angle_thlocal1)*pi/180 );
else
    angle_thlocal1 = angle_minthlocal1;
    cos_thlocal1 = cos_maxthlocal1;
end


% calculate error - difference between two iterations
squared_errors = sum((normal_preIter-faceNormals).^2, 2);
avg_error = sum(squared_errors) / numF;
error_largest = max( squared_errors );
fprintf('Iter: %d; largest error: %d; average error: %d\n', outerIter, error_largest, avg_error);



% record the previous-previous iteration
error_largest_pre = error_largest;
avg_error_pre = avg_error;


outerIter = outerIter + 1;
end
toc;

%% write to file: face normals
fid = fopen('FaceNormals.txt', 'w');
fprintf(fid, '%0.6f %0.6f %0.6f\n', [ faceNormals(:,1)'; faceNormals(:,2)'; faceNormals(:,3)' ]);
fclose(fid);


%% vertex update [Sun et al. 2007 TVCG paper]
tic;
newpos = zeros(numV, 3);
numberofIters = 10;% set by users
count_it = 0;
while (count_it < numberofIters)
    centroids = computeCentroidofTriangles(V, F);
    for i=1:numV
        sum = zeros(1,3);
        neF = vNf{i,1};
        for j=1:numel(neF)
            currF = neF(j);
            sum = sum + faceNormals(currF,:) * ( faceNormals(currF,:) * (centroids(currF,:)-V(i,:))' );
        end
        if numel(neF) > 0
            newpos(i,:) = V(i,:) + sum / numel(neF);
        end 
    end
    
    V = newpos;
    
    count_it = count_it + 1;
end
toc;
%write mesh ply file
Normv = zeros(numV,3);
writemeshPLY(V, F-1, Normv, 'denoisemesh.ply');




