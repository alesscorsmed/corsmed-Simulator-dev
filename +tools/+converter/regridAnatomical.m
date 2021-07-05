
clear all; close all; clc;

%% select model
modelName = 'XCAT';
switch lower(modelName)
    
    case 'xcat'
        modelPath  = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\originals\XCAT_model_20190326_234501.mat';
        pdPath     = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\originals\XCAT_PDinhom_20200506.mat';
        targetPath = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\new\XCAT_model_20201230';
        backgroundTissue = 79;
    case 'voxelman'
        modelPath  = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\originals\VoxelMan_model_20200520d.mat';
        pdPath     = [];
        targetPath = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\new\VoxelMan_model_20201230';
        backgroundTissue = 204;
    otherwise
        modelPath  = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\originals\MIDA_model_20200322.mat';
        pdPath     = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\originals\MIDA_PDinhom_20200506.mat';
        targetPath = 'C:\Users\jorge\CODE\EDUTOOL\modelGridder\new\MIDA_model_20201230';
        backgroundTissue = 50;
        
end

%% load model
model = load(modelPath);
if ~isempty(pdPath)
    pd = load(pdPath, 'PDinhom');
    hasPD = 1;
else
    pd = [];
    hasPD = 0;
end

%% find limits
scatterR        = model.model_spatial;
scatterVoxels   = size(scatterR,1);
% scatterPD       = model.inhom;

domain          = model.model_dim;
resolution      = model.model_elementsize;
tissueType      = model.model_new;
tissueValues    = model.model_tissues;
numTissues      = size(tissueValues,1);
numProperties   = size(tissueValues,2);

clear model;

minLimits = min(scatterR);
maxLimits = max(scatterR);

% linear discretization
x = minLimits(1):resolution(1):maxLimits(1);
y = minLimits(2):resolution(2):maxLimits(2);
z = minLimits(3):resolution(3):maxLimits(3);

% domain size
nX = length(x);
nY = length(y);
nZ = length(z);

% 3D domain
[xq,yq,zq] = meshgrid(x,y,z);

% sort scatterR in ascending z order
[~,idxSorted]   = sort(scatterR(:,3),'ascend');
scatterR(:,:)   = scatterR(idxSorted,:);
pd(:,1)         = pd(idxSorted,1);
tissueType(:,1) = tissueType(idxSorted,1);

clear idxSorted;

% loop on the z slices and perform scattered interpolation
for zz = 1:nZ

    tSlice = tic();
    
    % get corresponding query slice
    zSlice = unique(zq(:,:,zz));
    rqSlice = [ reshape(xq(:,:,zz),[],1), ...
        reshape(yq(:,:,zz),[],1), ...
        reshape(zq(:,:,zz),[],1) ];
    
    % restrict interpolation positions
    idxSlice = find( ( scatterR(:,3) >= zSlice-resolution(3) ) ...
             & ( scatterR(:,3) <= zSlice+resolution(3) ) );
         
    if isempty(idxSlice)
        vq = backgroundTissue*ones(nX*nY,1);
        pq = zeros(nX*nY,1);
    else
        if ( abs(unique(scatterR(idxSlice,3)) - zSlice) < 1e-8*zSlice )
            % 2D
            % create interpolant
            F = scatteredInterpolant( scatterR(idxSlice,1),...
                scatterR(idxSlice,2),...
                double(tissueType(idxSlice,1)), ...
                'nearest' );
            % evaluate interpolant
            vq = F( rqSlice(:,1), rqSlice(:,2) );
            
            if hasPD
                % create interpolant
                F = scatteredInterpolant( scatterR(idxSlice,1),...
                    scatterR(idxSlice,2),...
                    double(pd(idxSlice,1)), ...
                    'nearest' );
                % evaluate interpolant
                pq = F( rqSlice(:,1), rqSlice(:,2) );
            end
            
        else
            % 3D
            % create interpolant
            F = scatteredInterpolant( scatterR(idxSlice,1),...
                scatterR(idxSlice,2),...
                scatterR(idxSlice,3), ...
                double(tissueType(idxSlice,1)), ...
                'nearest' );
            % evaluate interpolant
            vq = F( rqSlice(:,1), rqSlice(:,2), rqSlice(:,3) );
            
            if hasPD
                % create interpolant
                F = scatteredInterpolant( scatterR(idxSlice,1),...
                    scatterR(idxSlice,2),...
                    scatterR(idxSlice,3), ...
                    double(pd(idxSlice,1)), ...
                    'nearest' );
                % evaluate interpolant
                pq = F( rqSlice(:,1), rqSlice(:,2), rqSlice(:,3) );
            end
            
        end
    end
    
    % reshape
    vq = reshape(vq,nY,nX,1);
    pq = reshape(pq,nY,nX,1);
    
    % remove nans
    vq(isnan(vq)) = backgroundTissue;
    pq(isnan(pq)) = 0.0;
    
    % save per slice
    save(sprintf('%s_%d.mat',targetPath,zz), 'vq', '-v7.3');
    if hasPD
        save(sprintf('%s_pd_%d.mat',targetPath,zz), 'pq', '-v7.3');
    end
    
    fprintf(1, '\n Slice %d interpolated: %.2f s', zz, toc(tSlice));
    
end
fprintf(1, '\n');
clear tissueType; clear pd;

spatial = [xq(:), yq(:), zq(:)];
clear xq; clear yq; clear zq;

%% create the new model
anatomicalModel = [];
anatomicalModel.name                = modelName;

% geometric characteristics
anatomicalModel.isGridded           = 1; % 1 if data is in meshgrid format
anatomicalModel.numVoxels           = nX*nY*nZ; % number of voxels
anatomicalModel.resolution          = resolution; % size of voxels [dx, dy, dz]
anatomicalModel.dimensions          = [nY, nX, nZ]; % grid dimensions: note matlab meshgrid permutes x and y
anatomicalModel.domain              = [ nX*resolution(1), nY*resolution(2), nZ*resolution(3) ]; % total domain size

% general properties
anatomicalModel.b0                  = 1.0;  % main field    : use to generate the initial magentization
anatomicalModel.mu                  = 1.0;  % susceptibility: use to generate the initial magnetization
anatomicalModel.pdFluctuationFactor = 0.05; % variation of PD in each tissue

% tisssues 
anatomicalModel.numTissues          = numTissues;
anatomicalModel.numProperties       = numProperties;
anatomicalModel.backgroundTissue    = backgroundTissue; % number of air tissue
anatomicalModel.fatTissuesIDs       = [ ]; % numbers of fat tissues
anatomicalModel.tissueValues        = tissueValues; % numTissues x numProperties
anatomicalModel.tissueDiff          = []; % Diffusion: numTissues x 3

% voxel values
anatomicalModel.spatial             = spatial; % numVoxels x 3 with coordinates: [x, y, z]
anatomicalModel.tissueType          = []; % numVoxels x 1 with tissue number
anatomicalModel.pdInhomogeneity     = []; % numVoxels x 1 with PD values
anatomicalModel.b0Inhomogeneity     = []; % numVoxels x 1 with B0 deltas


%% fill diffusion: values in the ballpark based on T1 and T2
dummyDiffusion = 1e-9*( 0.1 ...
    + 0.2*anatomicalModel.tissueValues(:,1) ...
    + 2*anatomicalModel.tissueValues(:,2) );
anatomicalModel.tissueDiff = zeros(anatomicalModel.numTissues,3);
anatomicalModel.tissueDiff(:,1) = 0.9*dummyDiffusion;
anatomicalModel.tissueDiff(:,2) = 1.0*dummyDiffusion;
anatomicalModel.tissueDiff(:,3) = 0.8*dummyDiffusion;


%% load and add properties
clear spatial;

%% tissue values
tissueType = zeros(nY,nX,nZ);
for zz = 1:nZ
    % save in case
    load(sprintf('%s_%d.mat',targetPath,zz));
    figure(1); imagesc(abs(vq)); colorbar;
    pause(0.01);
    % assign
    tissueType(:,:,zz) = vq;
    
end
anatomicalModel.tissueType = reshape(tissueType,[],1);
clear tissueType;

%% PD
if hasPD
    pdInhomogeneity = zeros(nY,nX,nZ);
    for zz = 1:nZ
        
        % save in case
        load(sprintf('%s_pd_%d.mat',targetPath,zz));
        figure(1); imagesc(abs(pq)); colorbar;
        pause(0.01);
        % assign
        pdInhomogeneity(:,:,zz) = pq;
        
    end
    anatomicalModel.pdInhomogeneity = reshape(pdInhomogeneity,[],1);
    clear pdInhomogeneity;
end


%% save the data
save(sprintf('%s.mat', targetPath), 'anatomicalModel', '-v7.3');

