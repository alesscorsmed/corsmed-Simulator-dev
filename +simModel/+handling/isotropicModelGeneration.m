function [sliceModel] = isotropicModelGeneration(slicePlane,...
    anatomicalModel, coilSystem, modelControl, dbgControl )
%
% SIMMODEL.HANDLING.ISOTROPICMODELGENERATION
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
functionName = 'simModel:handling:isotropicModelGeneration';
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

%% center in p1 (refPoint) and rotate the anatomical model
rotMat = (rotMatX*rotMatY*rotMatZ).';
if abs(p1) > 0
    refPoint = p1;
else
    refPoint = []; % reference point is origin
end

%% domain limits and discretization
% discretization step
dx = modelControl.gridStep(1);
dy = modelControl.gridStep(2);
dz = modelControl.gridStep(3);
% x domain discretization
xExtension = modelControl.xExtendPct; % extended percentage (def. 0.1)
xmin = -xExtension*p2Final(1); 
xmax = (1+xExtension)*p2Final(1);
x = union( -dx/2:-dx:xmin-dx, dx/2:dx:xmax+dx );
% y domain discretization
if strcmp(modelControl.coilMode,'basic')
    % Keep only the isochromats within the double FOV along the PE direction
    ymin = -0.5*p3Final(2); 
    ymax =  1.5*p3Final(2);
else % use largest possible dimensions (within certain limit)
    ymin = -2.0*p3Final(2); 
    ymax =  2.5*p3Final(2);
end
y = union( -dy/2:-dy:ymin-dy, dy/2:dy:ymax+dy );
% z domain discretization
if modelControl.useSliceThickness
    % take into account z thickness
    zmin = -(sliceThickness+dz)/2;
    zmax =  (sliceThickness+dz)/2;
    z = union( -dz/2:-dz:zmin, dz/2:dz:zmax );
else % force to 2D case    
    z = 0;
end

%% generate the meshgrid with target points
% query points on isotropic grid
[xq,yq,zq] = meshgrid(x,y,z);
rq = [ reshape(xq,[],1), reshape(yq,[],1), reshape(zq,[],1) ];

%% prepare the model, assign plane discretization data
sliceModel.r3D           = rq; % plane 3D grid in format: nx*ny*nz, 3
sliceModel.isotropic     = 1;
sliceModel.dimensions    = size(xq);  % dimensions of 3D plane grid (nx, ny, nz)
sliceModel.resolution    = [dx, dy, dz]; % voxel resolution (dx, dy, dz)

%% apply the interpolation for the anatomical model properties
if anatomicalModel.isGridded
    try
        %% anatomical model has a 3D structure that allows gridded interp
        [sliceModel] = simModel.interpolation.griddedAnatomical( ...
            sliceModel, rotMat, refPoint, anatomicalModel, dbgControl );
    catch
        if dbgControl.mode
            fprintf(fid, '\n%s : WARNING - gridded interpolation failed',...
                functionName);
        end
        %% if gridded interpolation fails, try scattered
        [sliceModel] = simModel.interpolation.scatteredAnatomical( ...
            sliceModel, rotMat, refPoint, anatomicalModel, dbgControl );
    end
else
    %% apply scattered interp
    [sliceModel] = simModel.interpolation.scatteredAnatomical( ...
        sliceModel, rotMat, refPoint, anatomicalModel, dbgControl );
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
    fprintf(fid, '\n%s : istropic model generated', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Isochromats     %d', sliceModel.numIsochromats);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end
