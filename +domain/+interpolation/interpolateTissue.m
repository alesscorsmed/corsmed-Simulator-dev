function [tissueMask] = interpolateTissue( rSlice, anatomicalModel )
%
% DOMAIN.INTERPOLATION.INTERPOLATETISSUE
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
functionName = 'domain.interpolation.interpolateTissue';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% prepare to do the interpolation in orginal 3D grid

% cropping: find limits of query points to restrict locations in 3D
rSliceMax = max(rSlice) + anatomicalModel.resolution;
rSliceMin = min(rSlice) - anatomicalModel.resolution;

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
    tissueMask = anatomicalModel.backgroundTissue*ones(size(rSlice,1),1);
    
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
    % nonZeroTissues = find( anatomicalModel.tissueValues(:,3)==0 );
    % idxZeroTissues = find( (isnan(vq)) & (ismember(vq,nonZeroTissues)) );
    idxZeroTissues = isnan(vq);
    vq(idxZeroTissues) = anatomicalModel.backgroundTissue;
    % assign to return
    tissueMask = vq;
    
end
