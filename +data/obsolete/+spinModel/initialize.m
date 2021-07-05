function [spinModel] = initialize(numSlices)
%
% DATA.SPINMODEL.INITIALIZE
%
%	Function that returns an empty
%   structure with a basic spin model.
%
%   This function is useful to define the fields.
%
% INPUT
%   numSlices   number of slices in the spinModel
%
% OUTPUT
%   spinModel   empty spinModel structure
%
%========================  CORSMED AB © 2020 ==============================
%

if (nargin < 1)
    numSlices = 1;
end
if isempty(numSlices) || (numSlices < 1)
    numSlices = 1;
end

% basic info
spinModel.name = 'EmptyModel';
spinModel.numSlices = numSlices;
spinModel.is3D = false;
spinModel.totalIsochromats = 0;

% spinModel.slice{1} has the first slice
% spinModel.slice{2} has the second slice
% ...
for ss = 1:spinModel.numSlices
    
    % initialize each slice plane coordinates and data
    %  (to be filled after plane handling)
    plane.LTop = 0.0; % Top-Left
    plane.RTop = 0.0; % Top-Right
    plane.LBot = 0.0; % Bottom-Left
    plane.RBot = 0.0; % Bottom-right
    % front end points
    plane.LTopFrontEnd  = 0.0; % Top-Left
    plane.RTopFrontEnd  = 0.0; % Top-Right
    plane.LBotFrontEnd  = 0.0; % Bottom-Left
    plane.RBotFrontEnd  = 0.0; % Bottom-right
    % transformations and orientations
    plane.p1OO    = 0.0; % center
    plane.p2FE    = 0.0; % Frequency Encoding axis
    plane.p3PE    = 0.0; % Phase Encoding axis
    plane.p2Final = 0.0; % Frequency Encoding axis final
    plane.p3Final = 0.0; % Phase Encoding axis final
    % rots and trans
    plane.rotMatX = zeros(3,3);
    plane.rotMatY = zeros(3,3);
    plane.rotMatZ = zeros(3,3);
    plane.transX = 0.0;
    plane.transY = 0.0;
    % coords and orientations
    plane.planeCoords = [];
    plane.BOrient = 0.0;
    plane.TOrient = 0.0;
    plane.LOrient = 0.0;
    plane.ROrient = 0.0;
    
    % initialize each slice model
    %  (to be filled after interpolation of anatomical model)
    model.r3D = []; % original 3D grid (in 3D grid format: nx, ny, nz, 3)
    model.resolution     = zeros(3,1); % voxel resolution (dx, dy, dz)
    
    model.nonZeroIndex   = []; % indexes with the entries of non zero
    model.numIsochromats = 0; % number of isochromats to simulate
    
    % positions of non-zero isochromats
    model.x  = zeros(model.numIsochromats,1); % vector with x positions
    model.y  = zeros(model.numIsochromats,1); % vector with y positions
    model.z  = zeros(model.numIsochromats,1); % vector with z positions
    
    % external voxel properties (not from anatomical model)
    model.bi         = zeros(model.numIsochromats,1);
    model.pd         = zeros(model.numIsochromats,1);
    model.xDiffusion = zeros(model.numIsochromats,1);
    model.yDiffusion = zeros(model.numIsochromats,1);
    model.zDiffusion = zeros(model.numIsochromats,1);
    
    % tissue types and property values
    model.numTissues    = 1; % number of tissue types for the model
    model.numProperties = 6; % number of tissue properties
    model.tissueValues  = zeros(model.numTissues,model.numProperties);
    model.tissueType    = zeros(model.numIsochromats,1);
    
    % coils
    model.numCoils  = 1;
    model.coilMapsX = zeros(model.numIsochromats,1);
    model.coilMapsY = zeros(model.numIsochromats,1);
    
    % update total number of useful isochromats in spinModel
    spinModel.totalIsochromats = spinModel.totalIsochromats ...
        + model.numIsochromats;
    
    % assign to slice and spinModel
    slice.plane = plane;
    slice.model = model;
    spinModel.slice{ss} = slice;

end

% initialize the slice for 3D cases
spinModel.slice3D = spinModel.slice{1};

