function TissueMaskResized = generateTissueMask(sliceModel,slicePlane,...
    acquisitionData,resizeFactorArray)

% % temporary solution to fix fail for anisotropic grids:
% % this avoids issues in other calls to this function 
% % besides from lin 426 of edt_run_RECONSTRUCTION.m
% if (nargin < 10) || isempty(struct_conf) || ~isfield(struct_conf,'isotropicGrid')
%     isotropicGrid = 1;
% else
%     isotropicGrid = struct_conf.isotropicGrid;
% end

uniqueZvalues = unique(sliceModel.z);
if size(uniqueZvalues,1)>1 %&& isotropicGrid % @@@
    
   % If there are more than 1 isochromats along the z direction 
   % of the anatomical model, keep only the one close to the xy
   % plane
   midSliceZ = uniqueZvalues(round(size(uniqueZvalues,1)/2),1);
   inds_mask = find(sliceModel.z==midSliceZ);
   
else
   inds_mask = 1:size(sliceModel.z,1);
   inds_mask = inds_mask';
end

anatModelIsotropicMaskStruct.x          = sliceModel.x(inds_mask,:);
anatModelIsotropicMaskStruct.y          = sliceModel.y(inds_mask,:);
anatModelIsotropicMaskStruct.z          = sliceModel.z(inds_mask,:);
anatModelIsotropicMaskStruct.tissueType = ...
   sliceModel.tissueType(inds_mask,:);

ind = find(...
    anatModelIsotropicMaskStruct.x(:,1)>=min(slicePlane.planeCoords(:,1))&...
    anatModelIsotropicMaskStruct.x(:,1)<=max(slicePlane.planeCoords(:,1))&...
    anatModelIsotropicMaskStruct.y(:,1)>=min(slicePlane.planeCoords(:,2))&...
    anatModelIsotropicMaskStruct.y(:,1)<=max(slicePlane.planeCoords(:,2)));

sliceInPlaneXYZ         = [anatModelIsotropicMaskStruct.x(ind,1),...
    anatModelIsotropicMaskStruct.y(ind,1),...
    anatModelIsotropicMaskStruct.z(ind,1)];
sliceInPlaneTissueType  = anatModelIsotropicMaskStruct.tissueType(ind);

% Bring the anatomical model back to original position where the p1 point
% is on (0,0,0), p2 point is on +x axis and p3 point is on +y axis.
sliceInPlaneXYZ(:,1) = sliceInPlaneXYZ(:,1) + slicePlane.transX;
sliceInPlaneXYZ(:,2) = sliceInPlaneXYZ(:,2) + slicePlane.transY;

% Create the isotropic grid
kspace(1,1) = acquisitionData.matrixX;
kspace(1,2) = acquisitionData.matrixY;

FOV(1,1) = acquisitionData.fovFE;
FOV(1,2) = acquisitionData.fovPE;

pixelSizeX = FOV(1,1)/kspace(1,1);
pixelSizeY = FOV(1,2)/kspace(1,2);

x = pixelSizeX/2:pixelSizeX:FOV(1,1);
y = pixelSizeY/2:pixelSizeY:FOV(1,2);
z = 0;

[X,Y,~] = meshgrid(x,y,z);

% @@@ Change this to interp3. It may not work with the currect version of
% XCAT
if(~isempty(sliceInPlaneXYZ))

    xqHover         = [X(:) Y(:)];
    sliceInPlane_py = transpose(sliceInPlaneXYZ(:,1:2));
    xqHover_py      = transpose(xqHover);
    F_tissue_IDs_py = py.scatteredNomatHover.scatteredPython(sliceInPlaneTissueType',sliceInPlane_py(:).',xqHover_py(:).',length(sliceInPlane_py),length(xqHover_py));
    F_tissue_IDs_py = double(py.array.array('d', py.numpy.nditer(F_tissue_IDs_py)));
    F_tissue_IDs    = F_tissue_IDs_py';

    tissuesIDs_newSlice     = F_tissue_IDs;
    tissuesIDs_newSlice_NEW = tissuesIDs_newSlice(:);
else
    tissue_size = kspace(1,2)*kspace(1,1);
    F_tissue_IDs = zeros(tissue_size,1);
    tissuesIDs_newSlice     = F_tissue_IDs;
    tissuesIDs_newSlice_NEW = tissuesIDs_newSlice(:);
end
    
TissueMask          = reshape(tissuesIDs_newSlice_NEW,kspace(1,2),kspace(1,1));
TissueMaskResized   = imresize(flipud(TissueMask),resizeFactorArray,'nearest');