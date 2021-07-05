function reconstructed_images = exportImage(contrastData,slicePlane,...
    expControl)

reconNx = size(contrastData.image,1);
reconNy = size(contrastData.image,1);

spatialCorners = [...
    num2str(slicePlane.LTopExtFOV(1)),',',num2str(slicePlane.LTopExtFOV(2)),',',num2str(slicePlane.LTopExtFOV(3)),';',...
    num2str(slicePlane.RTopExtFOV(1)),',',num2str(slicePlane.RTopExtFOV(2)),',',num2str(slicePlane.RTopExtFOV(3)),';',...
    num2str(slicePlane.LBotExtFOV(1)),',',num2str(slicePlane.LBotExtFOV(2)),',',num2str(slicePlane.LBotExtFOV(3)),';',...
    num2str(slicePlane.RBotExtFOV(1)),',',num2str(slicePlane.RBotExtFOV(2)),',',num2str(slicePlane.RBotExtFOV(3))];

%% Prepare text for image viewer
LTopText = '';
for i=1:length(contrastData.imageText.LTop)
    if ~strcmp(contrastData.imageText.LTop{i},'')
        LTopText = strcat(LTopText,'\n',contrastData.imageText.LTop{i});
    end
end

RTopText = '';
for i=1:length(contrastData.imageText.RTop)
    if ~strcmp(contrastData.imageText.RTop{i},'')
        RTopText = strcat(RTopText,'\n',contrastData.imageText.RTop{i});
    end
end

LBotText = '';
for i=1:length(contrastData.imageText.LBot)
    if ~strcmp(contrastData.imageText.LBot{i},'')
        LBotText = strcat(LBotText,'\n',contrastData.imageText.LBot{i});
    end
end

RBotText = '';
for i=1:length(contrastData.imageText.RBot)
    if ~strcmp(contrastData.imageText.RBot{i},'')
        RBotText = strcat(RBotText,'\n',contrastData.imageText.RBot{i});
    end
end

%% Add entry in db
if isfield(expControl,'approach') && ...
        strcmp(expControl.approach,'jsonstandalone')
    
    reconstructed_images.experiment_id      = num2str(expControl.experimentID);
    reconstructed_images.FOV_x              = num2str(slicePlane.fovX);
    reconstructed_images.FOV_y              = num2str(slicePlane.fovY);
    reconstructed_images.matrix_x           = num2str(reconNx);
    reconstructed_images.matrix_y           = num2str(reconNy);
    reconstructed_images.path               = strrep(contrastData.bmpName,'\','\\');
    reconstructed_images.spatial_corners    = spatialCorners;
    reconstructed_images.len_x              = num2str(slicePlane.lengthX);
    reconstructed_images.len_y              = num2str(slicePlane.lengthY);
    reconstructed_images.tl_info            = LTopText(3:end);
    reconstructed_images.tr_info            = RTopText(3:end);
    reconstructed_images.bl_info            = LBotText(3:end);
    reconstructed_images.br_info            = RBotText(3:end);
    reconstructed_images.download_path      = strrep(contrastData.dcmName,'\','\\');
    reconstructed_images.tissuemask_path    = strrep(contrastData.tissueMapBmpName,'\','\\');
    reconstructed_images.top_info           = slicePlane.TOrient;
    reconstructed_images.right_info         = slicePlane.ROrient;
    reconstructed_images.bottom_info        = slicePlane.BOrient;
    reconstructed_images.left_info          = slicePlane.LOrient;
    reconstructed_images.slice_seqnum       = num2str(contrastData.uniqueID);
    reconstructed_images.kspace_path        = strrep(contrastData.kspaceBmpName,'\','\\');
    
end