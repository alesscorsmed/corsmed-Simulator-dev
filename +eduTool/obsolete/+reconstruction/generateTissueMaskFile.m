function generateTissueMaskFile(outputImage,acquisitionData,sliceModel,...
    slicePlane,is3D,resizeFactorArray)

if sliceModel.numIsochromats~=0
    % If sliceModel.is3D = 1 pick that part of the slab of the anatomical 
    % model that corresponds to this slice. Exctract the tissue map from
    % this part
    if is3D
%         uniqueZvalues                       = unique(sliceModel.z);
%         noOfDiscreteSlices                  = size(uniqueZvalues,1);
%         noOfDiscreteSlicesPerAcqSlice       = noOfDiscreteSlices/totalSlicesWithinSlab;
%         noOfDiscreteSlicesTillFirstMiddle   = ceil(noOfDiscreteSlicesPerAcqSlice/2);
%         middleOfCurrentSlice = iSliceSlab*round(noOfDiscreteSlicesPerAcqSlice) + ...
%             noOfDiscreteSlicesTillFirstMiddle;
%         if middleOfCurrentSlice>noOfDiscreteSlices
%             middleOfCurrentSlice = noOfDiscreteSlices;
%         end
%         indsSliceMask = find(anatModelIsotropic_struct.model_spatial(:,3)==...
%             uniqueZvalues(middleOfCurrentSlice,1));
%         anatModelIsotropic4TissueMask_struct.model_spatial = ...
%             anatModelIsotropic_struct.model_spatial(indsSliceMask,:);
%         anatModelIsotropic4TissueMask_struct.coilmapsx     = ...
%             anatModelIsotropic_struct.coilmapsx(indsSliceMask,:);
%         anatModelIsotropic4TissueMask_struct.coilmapsy     = ...
%             anatModelIsotropic_struct.coilmapsy(indsSliceMask,:);
%         anatModelIsotropic4TissueMask_struct.model_new     = ...
%             anatModelIsotropic_struct.model_new(indsSliceMask,:);
    end

    tissueMask          = eduTool.reconstruction.generateTissueMask(...
        sliceModel,slicePlane,acquisitionData,resizeFactorArray);
    
    tissueMaskOutput    = domain.planeHandling.correctOrientation(...
        tissueMask,...
        slicePlane.LTopNew,slicePlane.RTopNew,...
        slicePlane.LBotNew,slicePlane.RBotNew,...
        slicePlane.TOrient,slicePlane.BOrient,...
        slicePlane.ROrient,slicePlane.LOrient,...
        acquisitionData.fovFE,acquisitionData.fovPE,...
        acquisitionData.foldoverDir);
    
    imwrite(uint8(tissueMaskOutput),outputImage.tissueMapBmpName);
end