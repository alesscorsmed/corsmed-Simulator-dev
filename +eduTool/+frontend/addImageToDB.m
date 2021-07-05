function addImageToDB(contrastData, slicePlane, sliceNum, expControl )

% addImageToDB(conn_localdb,experiment_id,thumbName,dicomName,...
%     structOfImage,sliceNumber,kspaceFilename,pulseq_id)
if isfield(expControl,'connLocalDB')
    sqlquery            = ['SELECT batch_id FROM pulse_sequence WHERE exper_id=',num2str(expControl.experimentID)];
    sqlquery_results    = exec(expControl.connLocalDB, sqlquery); 
    isBatch             = fetch(sqlquery_results);
else
    isBatch.Data{1,1} = NaN;
end


pulseq_id = expControl.pulseqID;
if(~isnan(isBatch.Data{1,1}) && exist('pulseq_id','var'))
    sqlquery        = ['SELECT file_id FROM batch_load WHERE id=',num2str(isBatch.Data{1,1})];
    sqlquery_results    = exec(expControl.connLocalDB, sqlquery);
    batchfile = fetch(sqlquery_results);
   
    sqlquery2 = ['SELECT pulseq_id,ROW_NUMBER() OVER (PARTITION BY batch_id ORDER BY id asc) AS row_num FROM batch_protocol WHERE batch_id =',num2str(isBatch.Data{1,1})];
    sql_results2 = exec(expControl.connLocalDB, sqlquery2);
    experiments = fetch(sql_results2);
    [rows,cols] = size(experiments.Data);
     
    for i=1:rows
        if expControl.pulseqID == experiments.Data{i,1}
        	experimentNum = experiments.Data{i,2};
        end
    end
        
    batchfile = extractAfter(batchfile.Data{1,1},'/Batch_Excel_');
    batchfile = extractBefore(batchfile,'.xlsx');

    destinationFolder = ['/efs-mount-point/Batches/Results/Batch_',...
        batchfile,'/',datestr(now,'dd_mm_yyyy'),'/Experiment_',num2str(experimentNum)];
    statusFolder = mkdir(destinationFolder); 
    statusCopy = copyfile(contrastData.bmpName,destinationFolder);
 
    destinationFolder2 = ['/efs-mount-point/Batches/Results/Batch_',...
        batchfile,'/CompareImages/Experiment_',num2str(experimentNum),'/Slice_',num2str(sliceNum)];
    statusFolder2 = mkdir(destinationFolder2); 
    statusCopy2 = copyfile(contrastData.bmpName,destinationFolder2);
    
end

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
if isfield(expControl,'connLocalDB')
    curs6 = exec(expControl.connLocalDB,...
        ['INSERT INTO edt_tool_local.reconstructed_images(',...
        'exper_id,pos_x,pos_y,pos_z,rot_x,rot_y,rot_z,FOV_x,',...
        'FOV_y,matrix_x,matrix_y,path,spatial_corners,len_x,len_y,course_id,',...
        'tl_info,tr_info,bl_info,br_info,download_path,tissuemask_path,',...
        'top_info,right_info,bottom_info,left_info,slice_seqnum,kspace_path) VALUES (',...
        num2str(expControl.experimentID),...
        ',',num2str(0),',',num2str(0),...
        ',',num2str(0),',',num2str(0),...
        ',',num2str(0),',',num2str(0),...
        ',',num2str(slicePlane.fovX),',',num2str(slicePlane.fovY),...
        ',',num2str(reconNx),',',num2str(reconNy),...
        ',''',strrep(contrastData.bmpName,'\','\\'),''',''',spatialCorners,''',',...
        num2str(slicePlane.lengthX),',',num2str(slicePlane.lengthY),',',...
        num2str(expControl.courseID),',''',...
        LTopText(3:end),''',''',RTopText(3:end),''',''',...
        LBotText(3:end),''',''',RBotText(3:end),''',''',...
        strrep(contrastData.dcmName,'\','\\'),''',''',...
        strrep(contrastData.tissueMapBmpName,'\','\\'),''',''',...
        slicePlane.TOrient,''',''',slicePlane.ROrient,''',''',...
        slicePlane.BOrient,''',''',slicePlane.LOrient,''',''',...
        num2str(contrastData.uniqueID),''',''',...
        strrep(contrastData.kspaceBmpName,'\','\\'),''')']);
elseif isfield(expControl,'applicationType') && ...
        strcmp(expControl.applicationType,'apijson')
    
end