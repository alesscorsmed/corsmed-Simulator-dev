function [expInfo] = experimentExecution(...
    anatomicalModel, coilSystem, mrSystem, ...
    acquisition, expControl, tagsStruct, parpoolConn, sessionData)
%
% EDUTOOL.RUN.EXPERIMENTEXECUTION
%
%	Runs a experiment.
%
% INPUT
%   sessionData         solution struct with initial data
%   anatomicalModel     struct with the anatomical model
%   coilModel           struct with models of the different coils
%
% OUTPUT
%   expInfo             experiment info for the user
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.experimentExecution';
if (nargin < 5)
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

%% check if Parpool is still on
p = eduTool.setup.openParpool();

%% try running simulation
try
    
    %% change status to 'started'
    eduTool.frontend.updateExperimentStatus(expControl.connLocalDB,expControl.experimentID,'started')
    
    %% create experiment progress columns and update to 1%
    exec(expControl.connLocalDB, ...
        ['INSERT INTO edt_tool_local.experiments_progress (exper_id,progress)',...
        ' VALUES (',num2str(expControl.experimentID),',1)']);

    %% SLICER   
    [simJobs, spinModel, pulseSequence, motionModel, sarReport, ...
    acquisition, expControl] = eduTool.run.slicer(...
    anatomicalModel, coilSystem, mrSystem, acquisition, expControl );
    
    %% check status
    eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)
    
    %% SIMULATOR
    [simJobs, expControl] = eduTool.run.engine( simJobs, ...
    pulseSequence, motionModel, expControl, tagsStruct, parpoolConn );
    
    %% check status
    eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)
    
    %% RECONSTRUCTION
    [reconData, expControl] = eduTool.run.recon(simJobs, ...
        spinModel, coilSystem, acquisition, expControl, sessionData,...
        pulseSequence, mrSystem);
    
    %% check status
    eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)
       
    %% DICON GENERATION
    %     %% process the image and send data back to FE
    %     [imageData] = dicom.processReconstruction( ...
    %         reconData, acquisition );
    
%     % plot for now plot
%     figure();
%     if spinModel.is3D
%         numZ = size(reconData.slice3D.iSpace,3);
%     else
%         numZ = spinModel.numSlices;
%     end
%     numC = reconData.numC;
%     count = 0;
%     for ci=1:numC
%         for zz=1:numZ
%             if spinModel.is3D
%                 iSpace = reconData.slice3D.iSpace;
%                 iMap = abs(iSpace(:,:,zz,1,ci).');
%             else
%                 iSpace = reconData.slice{zz}.iSpace;
%                 iMap = abs(iSpace(:,:,1,1,ci).');
%             end
%             % original image
%             count = count+1;
%             subplot(numC,numZ,count);
%             imagesc(iMap);
%             axis image; colormap gray; colorbar;
%             title(sprintf('slice %d contrast %d', zz, ci));
%         end
%     end
    pause(0.5); % so that it shows image

    %% update progress bar to full
    expControl.progress = 100;
    eduTool.frontend.progressUpdate(expControl);
    
    %% update status to 'finished'
    eduTool.frontend.updateExperimentStatus(expControl.connLocalDB,expControl.experimentID,'finished')

catch ME
    
    %% catch the error
    ME.identifier;
    ME.message;

    if strcmp(ME.message,'CANCELLED-BY-USER')
        %% if user cancellation
        formatOut = 'yyyy-mm-dd HH:MM:SS';
        timestamp = datestr(now,formatOut);
        errorMessage = sprintf(['%s - User (%d) - Exper (%d) - ',...
                ' Error Message:\n%s'], ...
                timestamp,expControl.userID,expControl.experimentID,ME.message);
    else               
        %% update status to 'Error'
        eduTool.frontend.updateExperimentStatus(expControl.connLocalDB,expControl.experimentID,'error')        
    end
    
    %% Write error message in db.
    % Do not write this message if this
    % comes from CX, JV and GB userID
    if ismember(expControl.userID,[790,933,1139])
        %% print error
        errorMessage = sprintf(['%s - User (%d) - Exper (%d) - ',...
            'Error in function %s() at line %d.',...
            '\n\nError Message:\n%s'], ...
            expControl.timeStamp, ...
            expControl.userID,...
            expControl.experimentID,...
            ME.stack(1).name,ME.stack(1).line,ME.message );
        fprintf(1, '%s\n', errorMessage);
        
    else % no developer
        
        %% save data for debugging
        errorDataFile = sprintf('%serrorBackend-%s-usr%s-exp%s.mat',...
            expControl.folderSystem.errorFolder, expControl.timeStamp, ...
            expControl.userID, expControl.experimentID);
        save(errorDataFile, 'expControl', 'acquisition', 'ME', '-v7.3');
        
        %% send error in connection to DB for backend
        errorMessage = sprintf(['%s - User (%d) - Exper (%d) - ',...
            'Error in function %s() at line %d.',...
            '\n\nError Message:\n%s',...
            '\n\nData for error replication saved in %s'], ...
            expControl.timeStamp, ...
            expControl.userID,...
            expControl.experimentID,...
            ME.stack(1).name,ME.stack(1).line,ME.message,...
            errorDataFile);
        fprintf(1, '%s\n', errorMessage);
        
        %% notify in slack
        eduTool.frontend.notifyAdminForErrors(expControl.connLocalDB,...
            errorMessage,expControl.instanceID)
    end
    
    %     %% inform user
    %     msg = ['The experiment stopped and quit unexpectedly. ',...
    %         'This issue has been reported to CORSMED for further review. ',...
    %         'If you were so kind as to help us in this process, ',...
    %         'please send us feedback (top right tab) with details about the experiment. ',...
    %         'We apologize and will try to address this issue in the coming release.'];
    %     eduTool.frontend.errorAndDBentry(sessionData.connLocalDB, msg, ...
    %         'cancelled-error',sessionData.experimentID,sessionData.pulseqID);
    
end

expInfo = 1;

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for experiment %d, elapsed time %.3fs',...
        functionName, expControl.experimentID, toc(tTotal));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
