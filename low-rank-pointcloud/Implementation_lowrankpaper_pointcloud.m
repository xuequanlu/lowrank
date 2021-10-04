%% structure, read data, and downsample?
%Point = struct('pos',{}, 'normal',{}, 'similarPoints',{}, 'neighborPoints',{});
[Position, Normal] = readPLY2('data/cube+1.0noise+PCA(15).ply');
[numPoints, dimPoints] = size(Position);

color = repmat([112, 138, 144, 255], numPoints, 1);

%bounding box
minBox = [1,1,1]; maxBox = [-1,-1,-1];
for i=1:numPoints
    if(minBox(1)>Position(i,1)) minBox(1) = Position(i,1); end
    if(minBox(2)>Position(i,2)) minBox(2) = Position(i,2); end
    if(minBox(3)>Position(i,3)) minBox(3) = Position(i,3); end
    if(maxBox(1)<Position(i,1)) maxBox(1) = Position(i,1); end
    if(maxBox(2)<Position(i,2)) maxBox(2) = Position(i,2); end
    if(maxBox(3)<Position(i,3)) maxBox(3) = Position(i,3); end
end
h0 = 4.0*sqrt( norm(maxBox-minBox)/numPoints );
%sometimes qradii is small to find ball neighbors (noise is large)
% this parameter is flexible and should be tuned 0.5h0~3.0h0
qradii = 2.0*h0; % should be big as we are finding non-local similar patches


%use built-in kdtree
tree = KDTreeSearcher( Position );% build the tree

%can be tuned: angle_thnonlocal, angle_minthnonlocal
angle_thnonlocal = 30.0
angle_minthnonlocal = 15.0
reduce_thnonlocal = 1.0;
cos_thnonlocal = cos(angle_thnonlocal*pi/180.0);
cos_maxthnonlocal = cos(angle_minthnonlocal*pi/180.0);

%-------representative normals for patches should be more robust than individual point normals---------
%can be tuned: angle_thlocal1, angle_minthlocal1
angle_thlocal1 = 30.0
angle_minthlocal1 = 15.0
reduce_thlocal1 = 1.0;
cos_thlocal1 = cos(angle_thlocal1*pi/180.0);
cos_maxthlocal1 = cos(angle_minthlocal1*pi/180.0);

reduce_rate = 1.1

%can be tuned: angle_tensor
angle_tensor = 30.0
angle_thlocal = angle_tensor/180.0*pi; %for text+wall: 10.0

%can be tuned: Patchsize_th, usually fixed
Patchsize_th = 1e10;


knn = 150 %a large enough region, fixed; default: 150

% Knn
[neighborPoints, dis_knn] = knnsearch(tree, Position, 'K', knn);%fast: increasing distance order
localknn = 60 %default 60, flexible to tune
localNNpoints = neighborPoints(:, 2:localknn);%remove itself
sigma_dis2 = ( 2.0*dis_knn(:, localknn) ).^2;%each point has a different radius


%% (1) initialization for normals: with anisotropic Tensor voting (optional)
fprintf('computing spatialW...');
tic;
spatialW = computeSpatialW(localNNpoints, Position, numPoints, sigma_dis2); % only compute once
toc;
init_iters = 1;%typically 1 iteration, since it is just initialization
fprintf('computing TensorVotingNormal for initial...');
tic;
Normal = TensorVotingNormal( localNNpoints, Normal, numPoints, init_iters, angle_thlocal, spatialW );
toc;

%global allFactors
allFactors = importdata('decomposedInteger.mat');

localIso = cell(1, numPoints);
simPatches = cell(1, numPoints);
%% the number of overall iterations: automatically determined by error rate

max_outerIter = 6;%max number of iterations
outerIter = 0;

tic;
while (outerIter<max_outerIter) 
%% (2) find similar patches for each point
%local representative normal of each patch, not the normal of this point!!!
   init_iters = 1;
   localstru = TensorVotingNormal( localNNpoints, Normal, numPoints, init_iters, angle_thlocal, spatialW );

parfor i = 1:numPoints
    idxs = neighborPoints(i, :); %enough large 'window'
    numofNeigh = numel(idxs);
    tmp_localstru = queryLocalstru(localstru, idxs);
    
    % record the isotropic neighbors consistent with localstructure
    idxs_ball = localNNpoints(i, :);% use local knn
    numofballNeigh = numel(idxs_ball);
    tmp_normal = queryNormal(Normal, idxs_ball);
    inAngle =  dot( repmat(localstru(i,:),numofballNeigh,1), tmp_normal, 2);
    
    [inAngle,idss] = sort(inAngle,'descend'); %decreasing order for cos(), increasing order in angle
    inAng1 = inAngle(inAngle >= cos_thlocal1);
    localIso{i} = [i, idxs_ball(:, idss(1:numel(inAng1),:) ) ];%should include the current point itself

    interAngle = dot( repmat(localstru(i,:),numofNeigh,1), tmp_localstru, 2) ;
    [interAngle, idxx] = sort(interAngle, 'descend');%increasing order in angle
    interAngle1 = interAngle(interAngle >= cos_thnonlocal);%angle_thnonlocal
    simPatches{i} = idxs( :, idxx(1:numel(interAngle1),:) ); %increasing error order, include itself
end

%% (3) for each point, perform nuclear norm minimization (low rank matrix approximation)
sumNormals = zeros(numPoints,3);
sumWeights = zeros(numPoints,1);

iter_svd = 2 % 1 iter is enough for real scanned point models, 2 is for CAD-like models
alpha_coeff = 1.0 %coefficient for weighting function, default 1.0
%record the abnormal point ids
abpoints = sparse(numPoints,1);
ccc = 0;
parfor i=1:numPoints
    idxs = simPatches{i};
    nonlocal_patchsize = numel(idxs);
    % by controling the number of non-local similar patches of patch i
    if nonlocal_patchsize > Patchsize_th
        nonlocal_patchsize = Patchsize_th;
        ccc = ccc + 1;
    end
    % acquire the similar patches for the i-th point; no need to re-sort all normal errors here
    tmpLocalIso = queryLocalIso( localIso, idxs(1,1:nonlocal_patchsize) );
    vec_tmpoints = cell2mat( tmpLocalIso );
    
    num_actual = numel(vec_tmpoints);
    
    if num_actual <=5 && outerIter>=2
        %handle abnormal points
        for tt=1:10% search 10 nearest neighbors
            tmInd = localNNpoints(i, tt);
            if ~any(tmInd==vec_tmpoints)
                abpoints(i,1) = tmInd;
                break;
            end
        end
    end
    
    [rowNum_patch, colNum_patch, Patchsize] = decomposeInteger( num_actual*3, allFactors );
    vec_tmpoints = vec_tmpoints(:,1:Patchsize);

    matVn_patch = queryNormal(Normal, vec_tmpoints);
    num_actual = numel(vec_tmpoints);
    
    
    matVn_patchX = matVn_patch;
  
    
    matVn_patchX = reshape(matVn_patchX, rowNum_patch, colNum_patch);
    
    % denoise each patch
    whichIter = 0;
    while whichIter < iter_svd 
        [U,SigmaX,V] = svd(matVn_patchX, 'econ');
        
        for kk=1:3
            if kk==1
                Temp = diag(SigmaX);
            end
            W_VecX = alpha_coeff * exp(- (2.0*Temp./Temp(1)).^2 ); 
            SigmaX1 = soft(SigmaX, diag(W_VecX));
            Temp = diag(SigmaX1);
        end

         X =  U*SigmaX1*V';
        
        % standardization for normals
        patch_n = reshape(X, [], 3);% N * 3
        patchnorms = repmat( sqrt(sum(patch_n.^2, 2)), 1, 3 );
        patch_n = patch_n ./ patchnorms;
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
    tmpPatchWeights = zeros(numPoints, 1);
    tmpPatchWeights(vec_true,:) = tmpPatchWeights(vec_true,:) + reco_id;
    sumWeights = sumWeights + tmpPatchWeights;
    %
    tmpPatchNormal = zeros(numPoints, 3);
    tmpPatchNormal(vec_true,:) = tmpPatchNormal(vec_true,:) + reco_n;
    sumNormals = sumNormals + tmpPatchNormal;
    
end


normal_preIter = Normal;%record the original normals
% average and standardization
Normal = sumNormals ./ ( repmat(sumWeights,1,3)  );% average, no influences
Allnorms = repmat( sqrt(sum(Normal.^2, 2)), 1, 3 )  ;
Normal = Normal ./ Allnorms;


%deal with specific abnormal points, after computing new normals
nnz_abpoints = nnz(abpoints);
if nnz_abpoints >= 1
    [rowc, colc, valuec] = find(abpoints);
    tmp = Normal(valuec, :);
    Normal(rowc, :) = tmp;
end

%%
fprintf('Iter: %d, angle_thnonlocal: %d; angle_thlocal1: %d\n', outerIter, angle_thnonlocal, angle_thlocal1);
fprintf('Iter: %d, cos_thnonlocal: %d; cos_thlocal1: %d\n', outerIter, cos_thnonlocal, cos_thlocal1);

%use cos() instead of angles, to reduce computation
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
squared_errors = sum((normal_preIter-Normal).^2, 2);
avg_error = sum(squared_errors) / numPoints;
error_largest = max( squared_errors );
fprintf('Iter: %d; largest error: %d; average error: %d\n', outerIter, error_largest, avg_error);

% record the previous-previous iteration
error_largest_pre = error_largest;
avg_error_pre = avg_error;

outerIter = outerIter + 1;
end

toc;

%% write to PLY file: only normal updated
color1 = color;
succ = writePLY( Position, Normal, color1,   strcat('output_normal_', sprintf('%d',outerIter-1), '.ply') );

%% point update [Lu et al 2020, Low Rank Matrix Approximation for 3D Geometry Filtering]
tic;
newpos = zeros(numPoints, 3);
numberofIters = 15;% set by users
count_it = 0;
kN = 15;% set by users
fold = 3.0;%by default
neighborPoints = neighborPoints(:,1:kN);
while (count_it < numberofIters)
    for i=1:numPoints
        sum = zeros(1,3);
        neP = neighborPoints(i,:);
        for j=1:numel(neP)
            t = neP(j);%neighboring point
            sum = sum + Normal(t,:) * ( Normal(t,:) * (Position(t,:)-Position(i,:))' );
            sum = sum + Normal(i,:) * ( Normal(i,:) * (Position(t,:)-Position(i,:))' );
        end
        if numel(neP) > 0
            newpos(i,:) = Position(i,:) + sum / (fold*numel(neP));
        end 
    end
    
    Position = newpos;
    
    count_it = count_it + 1;
end
toc;
%point position updated
succ = writePLY( Position, Normal, color1,   strcat('filteredPointCloud', '.ply') );




