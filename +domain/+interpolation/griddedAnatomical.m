function [model] = griddedAnatomical( model, rotMat, refPoint, ...
    anatomicalModel, dbgControl )
%
% DOMAIN.INTERPOLATION.GRIDDEDANATOMICAL
%
%   Interpolates for the query points ( positions of the model )
%   using the original anatomicalModel grid and data.
%
%   Query points in the model plane are rotated back to original
%   coordinates of the anatomicalModel, where the interpolation happens.
%
% INPUT
%
%
%
% OUTPUT
%   model    updated model struct with interpolated data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'domain.interpolation.griddedAnatomical';
if (nargin < 5)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
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

%% undo the transformations in the query points
% to bring them back to original coordinates
rSlice = model.r3D;
rSlice = rSlice * rotMat.'; % transpose of a rotation matrix brings back
if ~isempty(refPoint)
    rSlice = rSlice + refPoint; % translate back to ref point
end

% cropping: find limits of query points to restrict locations in 3D
rSliceMax = max(rSlice) + 2*anatomicalModel.resolution;
rSliceMin = min(rSlice) - 2*anatomicalModel.resolution;

%% prepare to do the interpolation in orginal 3D grid
% reshape anatomicalModel data to 3D
rAnatomical = reshape( anatomicalModel.spatial, [anatomicalModel.dimensions, 3] );

% valid indexes according to the cropping
idxX = find( ( rAnatomical(1,:,1,1) >= rSliceMin(1) ) ...
    & ( rAnatomical(1,:,1,1) <= rSliceMax(1) ) );
idxY = find( ( rAnatomical(:,1,1,2) >= rSliceMin(2) ) ...
    & ( rAnatomical(:,1,1,2) <= rSliceMax(2) ) );
idxZ = find( ( rAnatomical(1,1,:,3) >= rSliceMin(3) ) ...
    & ( rAnatomical(1,1,:,3) <= rSliceMax(3) ) );

% make sure there are at least 2 gridded points in each direction
% otherwise gridded interpolant may fail
if length(idxX) == 1
    if idxX == 1
        idxX = union(idxX, idxX+1);
    else
        idxX = union(idxX-1,idxX);
    end
end
if length(idxY) == 1
    if idxY == 1
        idxY = union(idxY, idxY+1);
    else
        idxY = union(idxY-1,idxY);
    end
end
if length(idxZ) == 1
    if idxZ == 1
        idxZ = union(idxZ, idxZ+1);
    else
        idxZ = union(idxZ-1,idxZ);
    end
end

%% apply the interpolation

if isempty(idxX) || isempty(idxY) || isempty (idxZ)
    
    %% slice is Out-of-bounds: Empty model
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
    %% slice is not Out-of-bounds
    
    % select domain within bounds, cropping
    x3d = rAnatomical(idxY,idxX,idxZ,1);
    y3d = rAnatomical(idxY,idxX,idxZ,2);
    z3d = rAnatomical(idxY,idxX,idxZ,3);
    
    % interpolate the tissueType
    v3d = reshape( anatomicalModel.tissueType, anatomicalModel.dimensions );
    v3d = v3d(idxY,idxX,idxZ);
    vq = interp3(x3d, y3d, z3d, v3d, ...
        rSlice(:,1), rSlice(:,2), rSlice(:,3), 'nearest');

    % find the indexes of useful tissues (non-zero tissues) and remove nans
    zeroTissues = find( anatomicalModel.tissueValues(:,3)==0 );
    idxTissues	   = find( (~isnan(vq)) & (~ismember(vq,zeroTissues)) );
    nIso           = length(idxTissues);
    
    % assign to model
    model.nonZeroIndex      = idxTissues; % indexes with the entries of non zero
    model.numIsochromats    = nIso; % number of isochromats to simulate
    model.tissueType        = reshape(vq(idxTissues), [nIso,1]);
    
    % interpolate the pdInhomogeneity if exists
    if ~isempty( anatomicalModel.pdInhomogeneity )
        v3d = reshape( anatomicalModel.pdInhomogeneity, anatomicalModel.dimensions );
        v3d = v3d(idxY,idxX,idxZ);
        vq  = interp3(x3d, y3d, z3d, v3d, ...
            rSlice(:,1), rSlice(:,2), rSlice(:,3), 'nearest');
        vq( isnan(vq) ) = 1.0;
        model.pd = reshape(vq(idxTissues), [nIso,1]);
    else
        model.pd = ones(nIso,1);
    end
    
    % interpolate the b0Inhomogeneity if exists
    if ~isempty( anatomicalModel.b0Inhomogeneity )
        v3d = reshape( anatomicalModel.b0Inhomogeneity, anatomicalModel.dimensions );
        v3d = v3d(idxY,idxX,idxZ);
        vq  = interp3(x3d, y3d, z3d, v3d, ...
            rSlice(:,1), rSlice(:,2), rSlice(:,3), 'linear');
        vq( isnan(vq) ) = 0.0;
        model.bi = reshape(vq(idxTissues), [nIso,1]);
    else
        model.bi = zeros(nIso,1);
    end
   
end

%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : gridded interpolation done', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Isochromats     %d', model.numIsochromats);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end