function [model] = anisotropicModelGeneration(plane,...
    anatomicalModel, coilSystem, expControl )
%
% DOMAIN.MODELHANDLING.ANISOTROPICMODELGENERATION
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
functionName = 'domain.anisotropicModelGeneration';
if (nargin < 4)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n\n%s : start', functionName);
end

%%  extract data from the slice plane a pre-compute data
p1              = plane.p1OO;
p2Final         = plane.p2Final;
p3Final         = plane.p3Final;
rotMatX         = plane.rotMatX;
rotMatY         = plane.rotMatY;
rotMatZ         = plane.rotMatZ;
transX          = plane.transX;
transY          = plane.transY;
sliceThickness  = plane.thickness;

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
xExtension = expControl.model.xExtendPct; % extended percentage (def. 0.1)
xmin = -xExtension*p2Final(1); 
xmax = (1+xExtension)*p2Final(1);
% y limits
if strcmp(expControl.model.coilMode,'basic')
    % Keep only the isochromats within the double FOV along the PE direction
    ymin = -0.5*p3Final(2); 
    ymax =  1.5*p3Final(2);
else % use largest possible dimensions (within certain limit)
    ymin = -2.0*p3Final(2); 
    ymax =  2.5*p3Final(2);
end
% z limits
if expControl.model.useSliceThickness
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
idxSlice = ( rAnatomical(:,1) >= xmin ) ...
    & ( rAnatomical(:,1) <= xmax ) ...
    & ( rAnatomical(:,2) >= ymin ) ...
    & ( rAnatomical(:,2) <= ymax ) ...
    & ( rAnatomical(:,3) >= zmin ) ...
    & ( rAnatomical(:,3) <= zmax );

%% assign data
if isempty(idxSlice)

    %% slice is Out-of-bounds: Empty model
    model.r3D           = []; % plane 3D grid in format: nx*ny*nz, 3
    model.isotropic     = 0;
    model.dimensions    = [];  % dimensions of 3D plane grid (nx, ny, nz)
    model.resolution    = [dx, dy, dz]; % voxel resolution (dx, dy, dz)
    model.nonZeroIndex   = []; % indexes with the entries of non zero
    model.numIsochromats = 0; % number of isochromats to simulate
    % external voxel properties (not from anatomical model)
    model.bi         = zeros(model.numIsochromats,1);
    model.pd         = zeros(model.numIsochromats,1);
    model.xDiffusion = zeros(model.numIsochromats,1);
    model.yDiffusion = zeros(model.numIsochromats,1);
    model.zDiffusion = zeros(model.numIsochromats,1);
    % tissue types
    model.tissueType = zeros(model.numIsochromats,1);
    
else
    
    %% slice is not Out-of-bounds: fill with data
    model.r3D            = rAnatomical(idxSlice,:); % plane 3D grid in format: nx*ny*nz, 3
    model.isotropic      = 0;
    model.dimensions     = [];  % dimensions of 3D plane grid (nx, ny, nz)
    model.resolution     = [dx, dy, dz]; % voxel resolution (dx, dy, dz)
    
    % tissue type
    vq = anatomicalModel.tissueType(idxSlice);
    
    % find the indexes of useful tissues (non-zero tissues) and remove nans
    nonZeroTissues = find( anatomicalModel.tissueValues(:,3)==0 );
    idxTissues	   = find( (~isnan(vq)) & (~ismember(vq,nonZeroTissues)) );
    nIso           = length(idxTissues);
    
    % assign to model
    model.nonZeroIndex      = idxTissues; % indexes with the entries of non zero
    model.numIsochromats    = nIso; % number of isochromats to simulate
    model.tissueType        = reshape(vq(idxTissues), [nIso,1]);
    
    % assign the pdInhomogeneity if exists
    if ~isempty( anatomicalModel.pdInhomogeneity )
        vq = anatomicalModel.pdInhomogeneity(idxSlice);
        vq( isnan(vq) ) = 1.0;
        model.pd = reshape(vq(idxTissues), [nIso,1]);
    else
        model.pd = ones(nIso,1);
    end
    
    % assign the b0Inhomogeneity if exists
    if ~isempty( anatomicalModel.b0Inhomogeneity )
        vq = anatomicalModel.b0Inhomogeneity(idxSlice);
        vq( isnan(vq) ) = 0.0;
        model.bi = reshape(vq(idxTissues), [nIso,1]);
    else
        model.bi = zeros(nIso,1);
    end
    
end

%% interpolate coil maps
[model] = domain.interpolation.griddedCoilMaps( ...
            model, rotMat, refPoint, coilSystem, expControl );

%% translate and assign coordinates
model.r3D(:,1) = model.r3D(:,1) - transX;
model.r3D(:,2) = model.r3D(:,2) - transY;
model.x = reshape(model.r3D(model.nonZeroIndex,1), [model.numIsochromats,1]);
model.y = reshape(model.r3D(model.nonZeroIndex,2), [model.numIsochromats,1]);
model.z = reshape(model.r3D(model.nonZeroIndex,3), [model.numIsochromats,1]);

%% assign not-interpolated data: tissue
model.b0            = anatomicalModel.b0;
model.mu            = anatomicalModel.mu;
model.numTissues    = anatomicalModel.numTissues; % number of tissue types for the model
model.numProperties = anatomicalModel.numProperties; % number of tissue properties
model.tissueValues  = anatomicalModel.tissueValues; % numTissues x numProperties with tissue values

%% rotate tissue diffussion values if exist
if ~isempty(anatomicalModel.tissueDiff)
    model.tissueDiff = anatomicalModel.tissueDiff*rotMat;
else
    model.tissueDiff = zeros(model.numTissues,3);
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : anistropic model generated', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Isochromats     %d', model.numIsochromats);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end
