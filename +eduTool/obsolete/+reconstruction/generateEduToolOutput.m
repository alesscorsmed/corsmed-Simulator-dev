function reconData = generateEduToolOutput(spinModel,acquisition,...
    sessionData,reconData,expControl,pulseSequence,mrSystem)

% Store images in the file system
if spinModel.is3D
    numZ = size(reconData.slice3D.iSpace,3);
else
    numZ = spinModel.numSlices;
end

% Initialize the dicom header
dicomStructInitial = eduTool.reconstruction.createDefaultDicomHeader();

numC    = reconData.numC;
count   = 0;
for ci=1:numC
    for zz=1:numZ
        if spinModel.is3D
            iSpace  = reconData.slice3D.iSpace;
            iMap    = abs(iSpace(:,:,zz,1,ci).');
            plane   = spinModel.slice3D.plane;
            model   = spinModel.slice3D.model;
        else
            iSpace  = reconData.slice{zz}.iSpace;
            iMap    = abs(iSpace(:,:,1,1,ci).');
            plane   = spinModel.slice{zz}.plane;
            model   = spinModel.slice{zz}.model;
        end
        
        sliceText = ['slice',num2str(zz),'_contrast',num2str(ci)];
        
        % Crop image along the FE so as to avoid foldover artifacts
        iMapCropped = ...
            reconstructor.image.cropImageAlongFE(iMap,...
            acquisition.reconstruction,acquisition.data.parallelImaging);
        
        imageUniqueID = (zz-1)*numZ + ci;
        
        % Resize images and save them in the file system
        reconData.slice{zz}.contrast{ci} = ...
            reconstructor.image.resizeImage(iMapCropped,...
            acquisition.reconstruction);
        
        % Create the folder where the results will be stored
        userFolderResults = [sessionData.folderSystem.gadgetronResultsFolder,...
            filesep,'user_',num2str(sessionData.userID)];
        if imageUniqueID == 1
            [s,cmdout] = system(['mkdir ',userFolderResults]);
        end

        % Create the folder where the reconstructed images will be stored
        outputFilenameNoExt = [sessionData.folderSystem.gadgetronResultsFolder,...
            filesep,'user_',num2str(sessionData.userID),...
            filesep,expControl.timestamp,'_',sliceText];

        % Change the orientation of the images. Calculate the 
        % XXXXOutputImage and save them in the structure reconData
        [outputImage,...
            LTopOutputImage,RTopOutputImage,...
            LBotOutputImage,RBotOutputImage,...
            topNew,bottomNew,rightNew,leftNew,...
            fovXNew,fovYNew,InPlanePhaseEncodingDirection] = ...
            domain.planeHandling.correctOrientation(...
            reconData.slice{zz}.contrast{ci}.normalizedImageResized,...
            plane.LTopNew,plane.RTopNew,plane.LBotNew,plane.RBotNew,...
            plane.TOrient,plane.BOrient,plane.ROrient,plane.LOrient,...
            acquisition.data.fovFE,acquisition.data.fovPE,...
            acquisition.data.foldoverDir);
        
        % Find the recon points for an extended FOV. The recon points for
        % the original image and .dcm files are already available in
        % plane.XXXXNew whereas the recon points for the image with the
        % correct orientation are stored in XXXXOutputImage.
        [LTopExtFOVOutputImage,RTopExtFOVOutputImage,...
            LBotExtFOVOutputImage,RBotExtFOVOutputImage] = ...
            eduTool.planeHandling.findReconPoints(LTopOutputImage,RTopOutputImage,...
            LBotOutputImage,RBotOutputImage,expControl.outerFOVratio);
        
        % Prepare text for image viewers
        textForImage = eduTool.frontend.prepareTextForImageViewer(...
            acquisition,pulseSequence,reconData,spinModel,expControl,...
            zz,ci);
        
        reconData.slice{zz}.contrast{ci}.outputImage.image          = outputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.LTop           = LTopOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.RTop           = RTopOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.LBot           = LBotOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.RBot           = RBotOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.LTopExtFOV     = LTopExtFOVOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.RTopExtFOV     = RTopExtFOVOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.LBotExtFOV     = LBotExtFOVOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.RBotExtFOV     = RBotExtFOVOutputImage;
        reconData.slice{zz}.contrast{ci}.outputImage.TOrient        = topNew;
        reconData.slice{zz}.contrast{ci}.outputImage.BOrient        = bottomNew;
        reconData.slice{zz}.contrast{ci}.outputImage.ROrient        = rightNew;
        reconData.slice{zz}.contrast{ci}.outputImage.LOrient        = leftNew;
        reconData.slice{zz}.contrast{ci}.outputImage.fovFE          = fovXNew;
        reconData.slice{zz}.contrast{ci}.outputImage.fovPE          = fovYNew;
        reconData.slice{zz}.contrast{ci}.outputImage.uniqueID       = imageUniqueID;
        reconData.slice{zz}.contrast{ci}.outputImage.imageText      = textForImage;
        reconData.slice{zz}.contrast{ci}.outputImage.foldoverDir    = InPlanePhaseEncodingDirection;

        tGenerateOutputTimer = tic;
        % Store the reconstructed images as .bmp and .txt files
        reconData.slice{zz}.contrast{ci}.outputImage.bmpName = ...
            sprintf('%s_thumbnail_corrOrient_MODULAR.bmp',outputFilenameNoExt);
        imwrite(reconData.slice{zz}.contrast{ci}.outputImage.image,...
            reconData.slice{zz}.contrast{ci}.outputImage.bmpName)
        reconData.slice{zz}.contrast{ci}.outputImage.txtName = ...
            sprintf('%s_thumbnail_corrOrient_MODULAR.txt',outputFilenameNoExt);
        tools.saveMatrixAsTxt(...
            reconData.slice{zz}.contrast{ci}.outputImage.image,...
            reconData.slice{zz}.contrast{ci}.outputImage.txtName);
        
        % Generate and store the .dcm file
        reconData.slice{zz}.contrast{ci}.outputImage.dcmName = ...
            sprintf('%s.dcm',outputFilenameNoExt);
        eduTool.reconstruction.generateDicomFile(...
            reconData.slice{zz}.contrast{ci}.outputImage,...
            acquisition,spinModel,expControl,...
            sessionData,mrSystem,dicomStructInitial)
        
        % Store the tissue mask
        % @@@ it does not work with 3D
        reconData.slice{zz}.contrast{ci}.outputImage.tissueMapBmpName = ...
            sprintf('%s_tissueMask_corrOrient.bmp',outputFilenameNoExt);
        eduTool.reconstruction.generateTissueMaskFile(...
            reconData.slice{zz}.contrast{ci}.outputImage,...
            acquisition.data,model,plane,spinModel.is3D,...
            reconData.slice{zz}.contrast{ci}.resizeFactorArray);
        
        % Store the kspace as .txt and .bmp files. This kspace is not the
        % initially created kspace. It is the kspace that comes from the FT
        % of the output image.
        reconData.slice{zz}.contrast{ci}.outputImage.kspaceTxtName = ...
            sprintf('%s_kspace.txt',outputFilenameNoExt);
        reconData.slice{zz}.contrast{ci}.outputImage.kspaceOutputComplex = ...
            reconData.operators.FT(...
            reconData.slice{zz}.contrast{ci}.outputImage.image);
        tools.saveMatrixAsTxt(...
            reconData.slice{zz}.contrast{ci}.outputImage.kspaceOutputComplex,...
            reconData.slice{zz}.contrast{ci}.outputImage.kspaceTxtName);
        reconData.slice{zz}.contrast{ci}.outputImage.kspaceBmpName = ...
            sprintf('%s_kspace.bmp',outputFilenameNoExt);
        reconData.slice{zz}.contrast{ci}.outputImage.kspaceOutputImage = ...
            mat2gray(log(abs(...
            reconData.slice{zz}.contrast{ci}.outputImage.kspaceOutputComplex)));
        imwrite(reconData.slice{zz}.contrast{ci}.outputImage.kspaceOutputImage,...
            reconData.slice{zz}.contrast{ci}.outputImage.kspaceBmpName);
        tGenerateOutputTimerDuration = toc(tGenerateOutputTimer);
        
        if expControl.debug.debugMode
            try % open file if possible, otherwise dump to stdout
                fid = fopen(expControl.debug.debugFile,'a');
            catch
                fid = 1;
            end
            
            fprintf(fid,...
                '\n  Elapsed Time for generation of output files %.3fs\n',...
                tGenerateOutputTimerDuration);
            
            if fid ~=1
                fclose(fid);
            end
        end
        
        % Add entry to db
        eduTool.frontend.addImageToDB(expControl,...
            reconData.slice{zz}.contrast{ci}.outputImage,zz)
    end
end