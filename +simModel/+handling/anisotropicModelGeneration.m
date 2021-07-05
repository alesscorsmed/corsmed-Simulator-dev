function [sliceModel] = anisotropicModelGeneration(slicePlane,...
    anatomicalModel, coilSystem, modelControl, dbgControl )
%
% SIMMODEL.HANDLING.ANISOTROPICMODELGENERATION
%
%     Function that generates the simulation model, by interpolating
%     data from the anatomical model into the disctretized
%     plane of the current slice
%
% INPUT
%
%
%
% OUTPUT
%   spinModel          updated spinModel struct with model
%
%========================  CORSMED AB © 2020 ==============================
%
%
functionName = 'simModel:handling:anisotropicModelGeneration';
% check args
if (nargin < 4) || isempty(slicePlane) || isempty(anatomicalModel) ...
        || isempty(coilSystem) || isempty(modelControl)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 5) || isempty(dbgControl)
    dbgControl.mode = 0;
    dbgControl.file = [];
end
% info for debugging
if dbgControl.mode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(dbgControl.file,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%%  extract data from the slice plane a pre-compute data
p1              = slicePlane.p1OO;
p2Final         = slicePlane.p2Final;
p3Final         = slicePlane.p3Final;
rotMatX         = slicePlane.rotMatX;
rotMatY         = slicePlane.rotMatY;
rotMatZ         = slicePlane.rotMatZ;
transX          = slicePlane.transX;
transY          = slicePlane.transY;
sliceThickness  = slicePlane.thickness;

% resolution step
dx = anatomicalModel.resolution(1);
dy = anatomicalModel.resolution(2);
dz = anatomicalModel.resolution(3);

%% center in p1 (refPoint) and rotate the anatomical model
rotMat = (rotMatX*rotMatY*rotMatZ).';
if abs(p1) > 0
    refPoint = p1;
else
    refPoint = []; % reference point is origin
end

%% find the limits to the domain
xExtension = modelControl.xExtendPct; % extended percentage (def. 0.1)
xmin = -xExtension*p2Final(1); 
xmax = (1+xExtension)*p2Final(1);
% y limits
if strcmp(modelControl.coilMode,'basic')
    % Keep only the isochromats within the double FOV along the PE direction
    ymin = -0.5*p3Final(2); 
    ymax =  1.5*p3Final(2);
else % use largest possible dimensions (within certain limit)
    ymin = -2.0*p3Final(2); 
    ymax =  2.5*p3Final(2);
end
% z limits
if modelControl.useSliceThickness
    % take into account z thickness
    zmin = -(sliceThickness + dz)/2;
    zmax =  (sliceThickness + dz)/2;
else % force to 2D case    
    zmin = -dz/2;
    zmax =  dz/2;
end

%% apply the transformations on the anatomical model points
% to bring them to the slice plane
rAnatomical = reshape( anatomicalModel.spatial, [], 3 );
if ~isempty(refPoint)
    rAnatomical = rAnatomical - refPoint; % translate to slice plane
end
rAnatomical = rAnatomical * rotMat; % rotate to the plane

% valid indexes according to the domain limits
idxSlice = find( ( rAnatomical(:,1) >= xmin ) ...
    & ( rAnatomical(:,1) <= xmax ) ...
    & ( rAnatomical(:,2) >= ymin ) ...
    & ( rAnatomical(:,2) <= ymax ) ...
    & ( rAnatomical(:,3) >= zmin ) ...
    & ( rAnatomical(:,3) <= zmax ) );

%% assign data
if isempty(idxSlice)

    %% slice is Out-of-bounds: Empty model
    sliceModel.r3D           = []; % plane 3D grid in format: nx*ny*nz, 3
    sliceModel.isotropic     = 0;
    sliceModel.dimensions    = [];  % dimensions of 3D plane grid (nx, ny, nz)
    sliceModel.resolution    = [dx, dy, dz]; % voxel resolution (dx, dy, dz)
    sliceModel.nonZeroIndex   = []; % indexes with the entries of non zero
    sliceModel.numIsochromats = 0; % number of isochromats to simulate
    % external voxel properties (not from anatomical model)
    sliceModel.bi         = zeros(sliceModel.numIsochromats,1);
    sliceModel.pd         = zeros(sliceModel.numIsochromats,1);
    sliceModel.xDiffusion = zeros(sliceModel.numIsochromats,1);
    sliceModel.yDiffusion = zeros(sliceModel.numIsochromats,1);
    sliceModel.zDiffusion = zeros(sliceModel.numIsochromats,1);
    % tissue types
    sliceModel.tissueType = zeros(sliceModel.numIsochromats,1);
    
else
    
    %% slice is not Out-of-bounds: fill with data
    sliceModel.r3D            = rAnatomical(idxSlice,:); % plane 3D grid in format: nx*ny*nz, 3
    sliceModel.isotropic      = 0;
    sliceModel.dimensions     = [];  % dimensions of 3D plane grid (nx, ny, nz)
    sliceModel.resolution     = [dx, dy, dz]; % voxel resolution (dx, dy, dz)

    % find the indexes of useful tissues: remove nans
    idxTissues = find(idxSlice(~isnan(anatomicalModel.tissueType(idxSlice)))); % pick the valid subset without NaNs
    nIso       = length(idxTissues);
    vq         = idxSlice(idxTissues); % relevant indexes of voxels to simulate
    
    % assign to model
    sliceModel.nonZeroIndex      = idxTissues; % indexes with the entries of non zero
    sliceModel.numIsochromats    = nIso; % number of isochromats to simulate
    
    if anatomicalModel.numFrames > 1
        % find the entries in the slice that are subject to frame motion
        [~,queryIdxFrame,origIdxFrame] = intersect(vq, anatomicalModel.idxFrame);
        % assign base tissue type
        vTissueType = reshape(anatomicalModel.tissueType(vq),nIso,1);
        % create a matrix of tissueTypes for each frame
        sliceModel.tissueType = repmat( vTissueType, [1, anatomicalModel.numFrames] );
        % modify the entries subject to motion
        sliceModel.tissueType(queryIdxFrame,:) = anatomicalModel.tissueTypeFrame(origIdxFrame,:);
    else
        % non motion tissue
        sliceModel.tissueType = reshape(anatomicalModel.tissueType(vq),nIso,1);
    end
    
    
    % assign the pdInhomogeneity if exists
    if ~isempty( anatomicalModel.pdInhomogeneity )
        sliceModel.pd = reshape(anatomicalModel.pdInhomogeneity(vq), [nIso,1]);
        sliceModel.pd(isnan(sliceModel.pd)) = 0.0;
    else
        % sliceModel.pd = ones(nIso,1);
        sliceModel.pd = []; % will later use the values from tissueValues
    end
    
    % assign the b0Inhomogeneity if exists
    if ~isempty( anatomicalModel.b0Inhomogeneity )
        sliceModel.bi = reshape(anatomicalModel.b0Inhomogeneity(vq), [nIso,1]);
        sliceModel.bi(isnan(sliceModel.bi)) = 0.0;
    else
        sliceModel.bi = zeros(nIso,1);
    end
    
end

%% interpolate coil maps
[sliceModel] = simModel.interpolation.griddedCoilMaps( ...
            sliceModel, rotMat, refPoint, coilSystem, dbgControl );

%% translate and assign coordinates
sliceModel.r3D(:,1) = sliceModel.r3D(:,1) - transX;
sliceModel.r3D(:,2) = sliceModel.r3D(:,2) - transY;
sliceModel.x = reshape(sliceModel.r3D(sliceModel.nonZeroIndex,1), [sliceModel.numIsochromats,1]);
sliceModel.y = reshape(sliceModel.r3D(sliceModel.nonZeroIndex,2), [sliceModel.numIsochromats,1]);
sliceModel.z = reshape(sliceModel.r3D(sliceModel.nonZeroIndex,3), [sliceModel.numIsochromats,1]);

%% assign not-interpolated data: tissue
sliceModel.b0            = anatomicalModel.b0;
sliceModel.mu            = anatomicalModel.mu;
sliceModel.numTissues    = anatomicalModel.numTissues; % number of tissue types for the model
sliceModel.numProperties = anatomicalModel.numProperties; % number of tissue properties
sliceModel.tissueValues  = anatomicalModel.tissueValues; % numTissues x numProperties with tissue values

%% rotate tissue diffussion values if exist
if ~isempty(anatomicalModel.tissueDiff)
    sliceModel.tissueDiff = anatomicalModel.tissueDiff*rotMat;
else
    sliceModel.tissueDiff = zeros(sliceModel.numTissues,3);
end

%% report
if dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : anistropic model generated', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Isochromats     %d', sliceModel.numIsochromats);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end
