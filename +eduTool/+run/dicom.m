function [imageData,jsonStr,jsonStructure] = dicom( ...
    reconData, imageData, expControl, sarReport)
%
% EDUTOOL.RUN.DICOM
%
%	Runs the reconstruction.
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.dicom';
if (nargin < 3)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% Initialize the dicom header
try
    dicomStructFile = 'inputs/dicomStructInitial.json';
    fjson = fopen(dicomStructFile,'r');
    dicomStructInitial = jsondecode(fread(fjson,inf,'*char').');
    fclose(fjson);
catch
    dicomStructInitial 	= image.dicom.createDefaultDicomHeader();
    dicomStructJson	= jsonencode(dicomStructInitial);
    dicomStructFile	= 'inputs/dicomStructInitial.json';
    fjson = fopen(dicomStructFile,'w');
    fwrite(fjson, dicomStructJson, 'char');
    fclose(fjson);
end

%% loop on the iSpace data and generate the images
% scaling factor: max abs in the stack of images
stackScale = 0;
% allocate space for multiple contrasts
multiContrastSlice = zeros(reconData.numX,reconData.numY,reconData.numC);
% loop
for sliceNum = 1:imageData.numSlices
    
    %% get the current plane
    slicePlane = imageData.slice{sliceNum}.plane;
    
    for frameNum = 1:imageData.numFrames
        
        %% get the contrasts
        if reconData.is3D
            multiContrastSlice(:,:,:) = reconData.slice{1}.frame{frameNum}.iSpace(:,:,sliceNum,:);
            tissueMask = reconData.slice{1}.mask(:,:,sliceNum);
        else
            multiContrastSlice(:,:,:) = reconData.slice{sliceNum}.frame{frameNum}.iSpace;
            tissueMask = reconData.slice{sliceNum}.mask;
        end
        
        %% process the contrasts to generate the images to store
        switch lower(imageData.processContrast)
            case 'waterfatoop'
                [cSpace, contrastText, isContrast]  = image.dicom.processWaterFatOOP( ...
                    multiContrastSlice, imageData);
            case 'swi'
                [cSpace, contrastText, isContrast]  = image.dicom.processSWI( ...
                    multiContrastSlice);
            case 'molli'
                [cSpace, contrastText, isContrast]  = image.dicom.processMOLLI( ...
                    multiContrastSlice, imageData.slice{sliceNum}.frame{frameNum}, ...
                    imageData.sequence.MOLLI, expControl);
            otherwise
                cSpace          = multiContrastSlice;
                contrastText    = [];
                isContrast      = ones(imageData.numContrasts,1);        
        end
        
        %% loop on the image contrasts (different from acquisition contrasts)
        for contrastNum = 1:imageData.numContrasts
            
            % extract the contrast data
            contrastData = imageData.slice{sliceNum}.frame{frameNum}.contrast{contrastNum};
            
            % update contrast image text, if required
            % line 3 is the one related to contrast
            textPosition = 3;
            if ~isempty(contrastText)
                contrastData.imageText.LTop{textPosition} = contrastText{contrastNum};
                textPosition = textPosition + 1;
            end
            % update frame text
            if imageData.numFrames > 1
                contrastData.imageText.LTop{textPosition} = ...
                    sprintf('Frame %d (CN %d, PN %d)', frameNum, ...
                    imageData.slice{sliceNum}.frame{frameNum}.ctIdx, ...
                    imageData.slice{sliceNum}.frame{frameNum}.phIdx);
            end
            
            
            % Resize and rotate images
            [contrastData] = image.dicom.processImage( contrastData, ...
                cSpace(:,:,contrastNum), tissueMask, slicePlane );
            
            % unique id
            contrastData.uniqueID = contrastNum ...
                + (sliceNum-1)*imageData.numContrasts*imageData.numFrames ...
                + (frameNum-1)*imageData.numContrasts;
            
            % define if it is a contrast (to scale) or a map
            if isContrast(contrastNum)
                contrastData.scaleImage = 1;
                % find max abs in image
                stackScale = max(stackScale, max(abs(contrastData.image(:))));
            else
                contrastData.scaleImage = 0;
            end
            
            %% assign to contrast
            imageData.slice{sliceNum}.frame{frameNum}.contrast{contrastNum} = contrastData;
            
        end
    end
end
       
% process and save the images
if stackScale > 0
    stackScale = 1/stackScale;
else
    stackScale = 0;
end

for sliceNum = 1:imageData.numSlices
    
    %% get the current plane
    slicePlane = imageData.slice{sliceNum}.plane;
    
    for frameNum = 1:imageData.numFrames
        
        %% loop on the image contrasts (different from acquisition contrasts)
        for contrastNum = 1:imageData.numContrasts
            
            %% save the data
            tGenerateOutputTimer = tic;
            
            %% extract contrast
            contrastData = imageData.slice{sliceNum}.frame{frameNum}.contrast{contrastNum};
            
            %% generate file names
            
            % slice text
            sliceText = ['slice',num2str(sliceNum), ...
                '_frame',num2str(frameNum), ...
                '_contrast',num2str(contrastNum)];
            
            % base file
            outputFilenameNoExt = sprintf('%s%s_%s',...
                expControl.folderSystem.experimentFolder, ...
                expControl.name, sliceText);
            
            % file type and image resolution
            if isfield(imageData, 'sequence') && isfield(imageData.sequence,'MOLLI')
                imageFileType   = '.png';
                imageResolution = 16;
                imageConverter  = @(image)uint16(image);
            else
                imageFileType = '.bmp';
                imageResolution = 8;
                imageConverter  = @(image)uint8(image);
            end
            
            % file names
            if isfield(expControl,'approach') && ...
                    strcmp(expControl.approach,'jsonstandalone')
                outputFilenameNoExtHighPriority = sprintf('%s%s_%s',...
                    expControl.folderSystem.highPriorityFolder, ...
                    expControl.name, sliceText);

                % Store the reconstructed normalized images as .bmp
                contrastData.bmpName = ...
                    sprintf('%s_thumbnail_corrOrient%s',...
                    outputFilenameNoExtHighPriority,imageFileType);               
            else
                % Store the reconstructed normalized images as .bmp
                contrastData.bmpName = ...
                    sprintf('%s_thumbnail_corrOrient%s',...
                    outputFilenameNoExt,imageFileType);
            end
            
            % store complex image as txt data
            contrastData.txtName = ...
                sprintf('%s_thumbnail_corrOrient.txt',outputFilenameNoExt);
            % Generate and store the .dcm file
            contrastData.dcmName = ...
                sprintf('%s.dcm',outputFilenameNoExt);
            % Store the tissue mask
            contrastData.tissueMapBmpName = ...
                sprintf('%s_tissueMask_corrOrient%s',outputFilenameNoExt,...
                imageFileType);
            % Store the kspace as .txt and .bmp files. This kspace is not the
            % initially created kspace. It is the kspace that comes from the FT
            % of the output image.
            contrastData.kspaceTxtName = ...
                sprintf('%s_kspace.txt',outputFilenameNoExt);
            % kSpace image
            contrastData.kspaceBmpName = ...
                sprintf('%s_kspace%s',outputFilenameNoExt,imageFileType);
            
            %% process data and save
            
            % scale image
            if contrastData.scaleImage
                % scale by stack max
                contrastData.imageScale = stackScale;
                contrastData.image      = stackScale*contrastData.image;
            else
                % self scale
                contrastData.imageScale = max(abs(contrastData.image(:)));
                if contrastData.imageScale > 0
                    contrastData.imageScale = 1/contrastData.imageScale;
                else
                    contrastData.imageScale = 0.0;
                end
                contrastData.image = contrastData.imageScale*contrastData.image;
            end
            % normalize: image will be already stack-scaled to [0,1] range
            contrastData.normalizedImage = abs(contrastData.image);
            % and convert to resolution
            contrastData.normalizedImage = imageConverter(...
                contrastData.normalizedImage*(2^imageResolution) );           
            
            %% write normalize image with required resolution
            if strcmpi(imageFileType, '.bmp')
                imwrite(contrastData.normalizedImage, ...
                    contrastData.bmpName );
            else % use bitdepth for resolution
                imwrite(contrastData.normalizedImage, ...
                    contrastData.bmpName, 'BitDepth', imageResolution);
            end
            
            %% write complex scaled image
            tools.saveMatrixAsTxt( contrastData.image, contrastData.txtName);
            
            %% generate the dicom
            % NOTE: we save the absolute of the image for the dicom
            image.dicom.generateDicomFile(...
                contrastData, slicePlane, imageData.bodyPartName, ...
                imageData.dicomInfo, dicomStructInitial);
            
            %% Store the tissue mask in required resolution
            % @@@ it does not work with 3D
            if strcmpi(imageFileType, '.bmp')
                imwrite(imageConverter(contrastData.mask), ...
                    contrastData.tissueMapBmpName );
            else % use bitdepth for resolution
                imwrite(imageConverter(contrastData.mask), ...
                    contrastData.tissueMapBmpName, 'BitDepth', imageResolution);
            end
            
            %% Store the kspace as .txt and .bmp files. 
            % This kspace is not the initially created kspace. 
            % It is the kspace that comes from the FT of the output image.
            contrastData.kspaceOutputComplex = ...
                reconData.encoding.operators.fFT( abs(contrastData.image) );
            % generate the image (log form) in the correct resolution
            contrastData.kspaceOutputImage = mat2gray( ...
                log(abs( contrastData.kspaceOutputComplex )) );
            contrastData.kspaceOutputImage = imageConverter( ...
                 contrastData.kspaceOutputImage *(2^imageResolution) );
            % save txt
            tools.saveMatrixAsTxt(...
                contrastData.kspaceOutputComplex, contrastData.kspaceTxtName);
            % save image of kspace
            if strcmpi(imageFileType, '.bmp')
                imwrite(contrastData.kspaceOutputImage, ...
                    contrastData.kspaceBmpName );
            else % use bitdepth for resolution
                 imwrite(contrastData.kspaceOutputImage, ...
                    contrastData.kspaceBmpName, 'BitDepth', imageResolution);
            end
            
            %% done
            tGenerateOutputTimerDuration = toc(tGenerateOutputTimer);
            
            if expControl.debug.debugMode
                fprintf(fid, '\n Files saved for slice %d contrast %d',...
                    sliceNum, contrastNum);
                fprintf(fid, '\n  Base file name  %s', outputFilenameNoExt);
                fprintf(fid, '\n  Elapsed Time    %.3fs', tGenerateOutputTimerDuration);
                fprintf(fid, '\n');
            end
            
            %% assign to contrast
            imageData.slice{sliceNum}.frame{frameNum}.contrast{contrastNum} = contrastData;
            
            % Add entry to db for eduTool app
            if strcmpi(expControl.application,'edutool')
                if isfield(expControl,'approach') && ...
                        strcmp(expControl.approach,'jsonstandalone')
                    
                    jsonStructure.reconstructed_images{contrastData.uniqueID} = ...
                        eduTool.frontend.exportImage(...
                        contrastData,slicePlane,expControl );
                    
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.tl_info = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.tl_info,'\n','/newline/');
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.tr_info = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.tr_info,'\n','/newline/');
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.bl_info = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.bl_info,'\n','/newline/');
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.br_info = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.br_info,'\n','/newline/');
                    
                    % Replace the path of the expControl.folderSystem.experimentFolder
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.path = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.path,...
                        expControl.folderSystem.highPriorityFolder,'');
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.download_path = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.download_path,...
                        expControl.folderSystem.experimentFolder,'');
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.tissuemask_path = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.tissuemask_path,...
                        expControl.folderSystem.experimentFolder,'');
                    jsonStructure.reconstructed_images{contrastData.uniqueID}.kspace_path = ...
                        strrep(jsonStructure.reconstructed_images{contrastData.uniqueID}.kspace_path,...
                        expControl.folderSystem.experimentFolder,'');
                    
                    jsonStr = jsonencode(jsonStructure);
                    
                else
                    
                    eduTool.frontend.addImageToDB( ...
                        contrastData, slicePlane, sliceNum, expControl );
                    
                    jsonStr = '';
                    
                end
            end
            
        end
    end
end

if strcmpi(expControl.application,'edutool') && ...
            isfield(expControl,'approach') && ...
            strcmp(expControl.approach,'jsonstandalone')
    % Add the path of extra files that accompany this experiment (eg. log)
    logDataFile = sprintf('%s%s_LOG.txt',expControl.folderSystem.logFolder, ...
        expControl.name);
    jsonStructure.log.log_path = ...
        strrep(logDataFile,expControl.folderSystem.experimentFolder,'');
    
    % Add the SAR and TIRL indicators
    jsonStructure.UIindicators = sarReport.UIindicators;
        
%     expControl.redis.R = tools.redis.redisSetJsonWrapper(expControl.redis.R,...
%         expControl.redis.keys.experimentResultsRedisKey,jsonStructure);
else
    jsonStructure = [];
end

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time          %.3fs', tTotal);
    fprintf(fid, '\n  Number of Slices      %d', imageData.numSlices);
    fprintf(fid, '\n  Number of Frames      %d', imageData.numFrames);
    fprintf(fid, '\n  Number of Contrasts   %d', imageData.numContrasts);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

% 
% 
% 
% 
% 
% function reconData = generateEduToolOutput(spinModel,acquisition,...
%     sessionData,reconData,expControl,pulseSequence,mrSystem)
% 
% % Store images in the file system
% if spinModel.is3D
%     numZ = size(reconData.slice3D.iSpace,3);
% else
%     numZ = spinModel.numSlices;
% end
% 
% % Initialize the dicom header
% dicomStructInitial = eduTool.reconstruction.createDefaultDicomHeader();
% 
% numC    = reconData.numC;
% count   = 0;
% for ci=1:numC
%     for zz=1:numZ
%         if spinModel.is3D
%             iSpace  = reconData.slice3D.iSpace;
%             iMap    = abs(iSpace(:,:,zz,1,ci).');
%             plane   = spinModel.slice3D.plane;
%             model   = spinModel.slice3D.model;
%         else
%             iSpace  = reconData.slice{zz}.iSpace;
%             iMap    = abs(iSpace(:,:,1,1,ci).');
%             plane   = spinModel.slice{zz}.plane;
%             model   = spinModel.slice{zz}.model;
%         end
%         
%         sliceText = ['slice',num2str(zz),'_contrast',num2str(ci)];
%         
%         % Crop image along the FE so as to avoid foldover artifacts
%         iMapCropped = ...
%             reconstructor.image.cropImageAlongFE(iMap,...
%             acquisition.reconstruction,acquisition.data.parallelImaging);
%         
%         imageUniqueID = (zz-1)*numZ + ci;
%         
%         % Resize images and save them in the file system
%         reconData.slice{zz}.contrast{ci} = ...
%             reconstructor.image.resizeImage(iMapCropped,...
%             acquisition.reconstruction);
%         
%         % Create the folder where the results will be stored
%         userFolderResults = [sessionData.folderSystem.gadgetronResultsFolder,...
%             filesep,'user_',num2str(sessionData.userID)];
%         if imageUniqueID == 1
%             [s,cmdout] = system(['mkdir ',userFolderResults]);
%         end
% 
%         % Create the folder where the reconstructed images will be stored
%         outputFilenameNoExt = [sessionData.folderSystem.gadgetronResultsFolder,...
%             filesep,'user_',num2str(sessionData.userID),...
%             filesep,expControl.timestamp,'_',sliceText];
% 
%         % Change the orientation of the images. Calculate the 
%         % XXXXOutputImage and save them in the structure reconData
%         [outputImage,...
%             LTopOutputImage,RTopOutputImage,...
%             LBotOutputImage,RBotOutputImage,...
%             topNew,bottomNew,rightNew,leftNew,...
%             fovXNew,fovYNew,InPlanePhaseEncodingDirection] = ...
%             domain.planeHandling.correctOrientation(...
%             reconData.slice{zz}.contrast{ci}.normalizedImageResized,...
%             plane.LTopNew,plane.RTopNew,plane.LBotNew,plane.RBotNew,...
%             plane.TOrient,plane.BOrient,plane.ROrient,plane.LOrient,...
%             acquisition.data.fovFE,acquisition.data.fovPE,...
%             acquisition.data.foldoverDir);
%         
%         % Find the recon points for an extended FOV. The recon points for
%         % the original image and .dcm files are already available in
%         % plane.XXXXNew whereas the recon points for the image with the
%         % correct orientation are stored in XXXXOutputImage.
%         [LTopExtFOVOutputImage,RTopExtFOVOutputImage,...
%             LBotExtFOVOutputImage,RBotExtFOVOutputImage] = ...
%             eduTool.planeHandling.findReconPoints(LTopOutputImage,RTopOutputImage,...
%             LBotOutputImage,RBotOutputImage,expControl.outerFOVratio);
%         
%         % Prepare text for image viewers
%         textForImage = eduTool.frontend.prepareTextForImageViewer(...
%             acquisition,pulseSequence,reconData,spinModel,expControl,...
%             zz,ci);
%         
%         imageData.slice{sliceNum}.contrast{contrastNum}.image          = outputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.LTop           = LTopOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.RTop           = RTopOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.LBot           = LBotOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.RBot           = RBotOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.LTopExtFOV     = LTopExtFOVOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.RTopExtFOV     = RTopExtFOVOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.LBotExtFOV     = LBotExtFOVOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.RBotExtFOV     = RBotExtFOVOutputImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.TOrient        = topNew;
%         imageData.slice{sliceNum}.contrast{contrastNum}.BOrient        = bottomNew;
%         imageData.slice{sliceNum}.contrast{contrastNum}.ROrient        = rightNew;
%         imageData.slice{sliceNum}.contrast{contrastNum}.LOrient        = leftNew;
%         imageData.slice{sliceNum}.contrast{contrastNum}.fovFE          = fovXNew;
%         imageData.slice{sliceNum}.contrast{contrastNum}.fovPE          = fovYNew;
%         imageData.slice{sliceNum}.contrast{contrastNum}.uniqueID       = imageUniqueID;
%         imageData.slice{sliceNum}.contrast{contrastNum}.imageText      = textForImage;
%         imageData.slice{sliceNum}.contrast{contrastNum}.foldoverDir    = InPlanePhaseEncodingDirection;
% 
%         tGenerateOutputTimer = tic;
%         % Store the reconstructed images as .bmp and .txt files
%         imageData.slice{sliceNum}.contrast{contrastNum}.bmpName = ...
%             sprintf('%s_thumbnail_corrOrient.bmp',outputFilenameNoExt);
%         imwrite(imageData.slice{sliceNum}.contrast{contrastNum}.image,...
%             imageData.slice{sliceNum}.contrast{contrastNum}.bmpName)
%         imageData.slice{sliceNum}.contrast{contrastNum}.txtName = ...
%             sprintf('%s_thumbnail_corrOrient.txt',outputFilenameNoExt);
%         tools.saveMatrixAsTxt(...
%             imageData.slice{sliceNum}.contrast{contrastNum}.image,...
%             imageData.slice{sliceNum}.contrast{contrastNum}.txtName);
%         
%         % Generate and store the .dcm file
%         imageData.slice{sliceNum}.contrast{contrastNum}.dcmName = ...
%             sprintf('%s.dcm',outputFilenameNoExt);
%         eduTool.reconstruction.generateDicomFile(...
%             reconData.slice{zz}.contrast{ci}.outputImage,...
%             acquisition,spinModel,expControl,...
%             sessionData,mrSystem,dicomStructInitial)
%         
%         % Store the tissue mask
%         % @@@ it does not work with 3D
%         imageData.slice{sliceNum}.contrast{contrastNum}.tissueMapBmpName = ...
%             sprintf('%s_tissueMask_corrOrient.bmp',outputFilenameNoExt);
%         eduTool.reconstruction.generateTissueMaskFile(...
%             reconData.slice{zz}.contrast{ci}.outputImage,...
%             acquisition.data,model,plane,spinModel.is3D,...
%             reconData.slice{zz}.contrast{ci}.resizeFactorArray);
%         
%         % Store the kspace as .txt and .bmp files. This kspace is not the
%         % initially created kspace. It is the kspace that comes from the FT
%         % of the output image.
%         imageData.slice{sliceNum}.contrast{contrastNum}.kspaceTxtName = ...
%             sprintf('%s_kspace.txt',outputFilenameNoExt);
%         imageData.slice{sliceNum}.contrast{contrastNum}.kspaceOutputComplex = ...
%             reconData.operators.FT(...
%             imageData.slice{sliceNum}.contrast{contrastNum}.image);
%         tools.saveMatrixAsTxt(...
%             imageData.slice{sliceNum}.contrast{contrastNum}.kspaceOutputComplex,...
%             imageData.slice{sliceNum}.contrast{contrastNum}.kspaceTxtName);
%         imageData.slice{sliceNum}.contrast{contrastNum}.kspaceBmpName = ...
%             sprintf('%s_kspace.bmp',outputFilenameNoExt);
%         imageData.slice{sliceNum}.contrast{contrastNum}.kspaceOutputImage = ...
%             mat2gray(log(abs(...
%             imageData.slice{sliceNum}.contrast{contrastNum}.kspaceOutputComplex)));
%         imwrite(imageData.slice{sliceNum}.contrast{contrastNum}.kspaceOutputImage,...
%             imageData.slice{sliceNum}.contrast{contrastNum}.kspaceBmpName);
%         tGenerateOutputTimerDuration = toc(tGenerateOutputTimer);
%         
%         if expControl.debug.debugMode
%             try % open file if possible, otherwise dump to stdout
%                 fid = fopen(expControl.debug.debugFile,'a');
%             catch
%                 fid = 1;
%             end
%             
%             fprintf(fid,...
%                 '\n  Elapsed Time for generation of output files %.3fs\n',...
%                 tGenerateOutputTimerDuration);
%             
%             if fid ~=1
%                 fclose(fid);
%             end
%         end
%         
%         % Add entry to db
%         eduTool.frontend.addImageToDB(expControl,...
%             reconData.slice{zz}.contrast{ci}.outputImage)
%     end
% end