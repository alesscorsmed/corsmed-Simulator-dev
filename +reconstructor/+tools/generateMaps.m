function [tissueMask,Cx,Cy] = generateMaps( plane, encodingPlan, ...
    anatomicalModel, coilSystem, expControl )
%
% RECONSTRUCTOR.TOOLS.GENERATEMAPS
%
%   Interpolates for the query points ( positions of the model )
%   using the original Coil grid and data.
%
%   Query (slice) points are rotated back to the coil coordinates.
%   The interpolation is applied on those original coordinates,
%   assuming a structured coordinate system (gridded interpolation).
%
%
% INPUT
%     model         struct with slice anatomical model data 
%     rotMat        rotation matrix
%     refPoint      reference point for rotation
%     coilSystem    struct with coils data
%     expControl    control info
%
% OUTPUT
%     model    updated model struct with interpolated data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructor.tools.generateMaps';
if (nargin < 5)
    ME = MException('Sense:wrongArgCount',...
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
    fprintf(fid, '\n%s : start', functionName);
end

%% generate query points from the slice
% Generate the slice grid in Image coordinates: query points
nX = encodingPlan.imSizeX;
nY = encodingPlan.imSizeY;
nZ = encodingPlan.imSizeZ;
% get FOV of image
fovX = encodingPlan.imFovX;
fovY = encodingPlan.imFovY;
% fovZ = encodingPlan.imFovZ;
% pixel size
pixelSizeX = fovX/nX;
pixelSizeY = fovY/nY;
% pixelSizeZ = fovZ/nZ;
% plane so that p1 = (0,0), p2 is on FE axis, p3 is on PE axis
x = pixelSizeX/2:pixelSizeX:fovX;
y = pixelSizeY/2:pixelSizeY:fovY;
if nZ > 1
    % coordinates as they come from the rotated planes
    z = encodingPlan.zPositions;
    % % z coordinates around 0
    % if mod(nZ,2) % odd number, include 0
    %     z = 0:pixelSizeZ:fovZ/2;
    % else % even number, symmetric around zero
    %     z = pixelSizeZ/2:pixelSizeZ:fovZ/2;
    % end
    % z = union(z, -z);
else
    % on plane if single slice
    z = 0;
end
% mesh grid of points
[xq,yq,zq] = meshgrid(x,y,z);

%% undo the domain transformations 
% in order to move for XY plane to original positions in 3D
% extract data from the slice domain
refPoint    = plane.p1OO;
rotMatX     = plane.rotMatX;
rotMatY     = plane.rotMatY;
rotMatZ     = plane.rotMatZ;
% rotMat is already transposed: gets us from plane to original
rotMat = (rotMatX*rotMatY*rotMatZ); 
% undo translation and rotation
rSlice = [ reshape(xq,[],1), reshape(yq,[],1), reshape(zq,[],1) ];
rSlice = rSlice * rotMat;
if ~isempty(refPoint)
    rSlice = rSlice + refPoint; % translate back to ref point
end

%% tissue Mask interpolation
[tissueMask] = domain.interpolation.interpolateTissue( ...
    rSlice, anatomicalModel );
% reshape
tissueMask = reshape(tissueMask,nY,nX,nZ);
tissueMask = permute(tissueMask,[2,1,3]);
tissueMask = flip(tissueMask,2);

%% parallel receive: Use the Bx and By maps to create the B1- maps
% get the active receive coil structure
coilStruct = coilSystem.coilModel{coilSystem.indexRx};
% interpolate maps into slice
rxSensType = 'pRX';
[Cx,Cy,rxSensType] = coils.interpolateCoilSens( ...
    rSlice, coilStruct, rxSensType, ...
    coilSystem.isocenter, coilSystem.b1mScaling );
% reshape into 3D
% Notice: matlab interp3 does y,x,z interpolation for Coils
% so reorganize to make sense
rxNumCoils = size(Cx,2);
Cx = reshape(Cx,nY,nX,nZ,rxNumCoils);
Cx = permute(Cx,[2,1,3,4]);
Cx = flip(Cx,2);
Cy = reshape(Cy,nY,nX,nZ,rxNumCoils);
Cy = permute(Cy,[2,1,3,4]);
Cy = flip(Cy,2);

%% report
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : mask and coil maps generated for %s',...
        functionName, coilStruct.data.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Query Points    %d', size(rSlice,1));
    fprintf(fid, '\n  RX Coil (#%2d)     %s', ...
        coilSystem.indexRx, ...
        coilSystem.coilModel{coilSystem.indexRx}.data.name);
    fprintf(fid, '\n  RX op. mode       %s', rxSensType);
    fprintf(fid, '\n  RX # Coil Maps    %d', rxNumCoils);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

