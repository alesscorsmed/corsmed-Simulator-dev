function [model] = isotropicModelGeneration(plane,...
    anatomicalModel, coilSystem, expControl )
%
% DOMAIN.GENERATESIMULATIONMODEL
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
functionName = 'domain.isotropicModelGeneration';
if (nargin < 1)
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

%% center in p1 (refPoint) and rotate the anatomical model
rotMat = (rotMatX*rotMatY*rotMatZ).';
if abs(p1) > 0
    refPoint = p1;
else
    refPoint = []; % reference point is origin
end

%% domain limits and discretization
% discretization step
dx = expControl.model.gridStep(1);
dy = expControl.model.gridStep(2);
dz = expControl.model.gridStep(3);
% x domain discretization
xExtension = expControl.model.xExtendPct; % extended percentage (def. 0.1)
xmin = -xExtension*p2Final(1); 
xmax = (1+xExtension)*p2Final(1);
x = union( -dx/2:-dx:xmin-dx, dx/2:dx:xmax+dx );
% y domain discretization
if strcmp(expControl.model.coilMode,'basic')
    % Keep only the isochromats within the double FOV along the PE direction
    ymin = -0.5*p3Final(2); 
    ymax =  1.5*p3Final(2);
else % use largest possible dimensions (within certain limit)
    ymin = -2.0*p3Final(2); 
    ymax =  2.5*p3Final(2);
end
y = union( -dy/2:-dy:ymin-dy, dy/2:dy:ymax+dy );
% z domain discretization
if expControl.model.useSliceThickness
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
model.r3D           = rq; % plane 3D grid in format: nx*ny*nz, 3
model.isotropic     = 1;
model.dimensions    = size(xq);  % dimensions of 3D plane grid (nx, ny, nz)
model.resolution    = [dx, dy, dz]; % voxel resolution (dx, dy, dz)

%% apply the interpolation for the anatomical model properties
if anatomicalModel.isGridded
    try
        %% anatomical model has a 3D structure that allows gridded interp
        [model] = domain.interpolation.griddedAnatomical( ...
            model, rotMat, refPoint, anatomicalModel, expControl );
    catch
        if expControl.debug.debugMode
            fprintf(fid, '\n%s : WARNING - gridded interpolation failed',...
                functionName);
        end
        %% if gridded interpolation fails, try scattered
        [model] = domain.interpolation.scatteredAnatomical( ...
            model, rotMat, refPoint, anatomicalModel, expControl );
    end
else
    %% apply scattered interp
    [model] = domain.interpolation.scatteredAnatomical( ...
        model, rotMat, refPoint, anatomicalModel, expControl );
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
    fprintf(fid, '\n%s : istropic model generated', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Isochromats     %d', model.numIsochromats);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end
