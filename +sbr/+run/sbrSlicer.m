function [spinModel,motionModel] = sbrSlicer(xMin,xMax,yMin,yMax,...
    zMin,zMax,gridStep,tissuesX,tissuesY,tissueValues)
%
% sbr.run.sbrSlicer
%
%	Creates the anatomical model with tissuesX x tissuesY tissues within
%	the single-slice anatomical model.
%
% INPUT
%   gridStep is a [1x3] array representing the size of the element in 3D
%
%
% OUTPUT
%

%% Create a grid, fill in the model_spatial
x = xMin:gridStep(1,1):xMax;
y = yMin:gridStep(1,2):yMax;
z = zMin:gridStep(1,3):zMax;

[X,Y,Z] = meshgrid(x,y,z);

model_spatial = [X(:),Y(:),Z(:)];

%% Fill in the model_new
% Label isochromats based on where they are located within a 
% tissuesX x tissuesY grid
model_new = ones(size(model_spatial,1),1);

x_breaks = linspace(xMin,xMax,(tissuesX+1));
y_breaks = linspace(yMin,yMax,(tissuesY+1));

tissueID = 0;
for i=1:length(x_breaks)-1
    tissueIndecesX = find(model_spatial(:,1)>=x_breaks(i) & ...
        model_spatial(:,1)<=x_breaks(i+1));
    for j=1:length(y_breaks)-1
        tissueID = tissueID + 1;
        tissueIndecesY = find(model_spatial(:,2)>=y_breaks(j) & ...
            model_spatial(:,2)<=y_breaks(j+1));
        tissueIndeces = intersect(tissueIndecesX,tissueIndecesY);
        model_new(tissueIndeces)=tissueID;
    end
end

%% Develop the spinModel
spinModel.numIsochromats    = size(model_spatial,1);
spinModel.resolution        = gridStep;
spinModel.mu                = 1;
spinModel.b0                = 1;
spinModel.x                 = model_spatial(:,1);
spinModel.y                 = model_spatial(:,2);
spinModel.z                 = model_spatial(:,3);
spinModel.bi                = zeros(size(model_spatial,1),1);
spinModel.pd                = ones(size(model_spatial,1),1);
spinModel.tissueDiff        = zeros(size(model_spatial,1),3);
spinModel.tissueType        = model_new;
spinModel.rxCoilMapsX       = ones(size(model_spatial,1),1);
spinModel.rxCoilMapsY       = zeros(size(model_spatial,1),1);
spinModel.numRxCoils        = 1;

if nargin<10
    spinModel.tissueValues  = zeros(tissueID,6);
else
    spinModel.tissueValues  = tissueValues;
end

%% Motion model
motionModel.type            = 'none';