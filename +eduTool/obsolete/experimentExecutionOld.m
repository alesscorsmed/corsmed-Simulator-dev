function [expInfo] = experimentExecution(...
    anatomicalModel, coilSystem, mrSystem, ...
    acquisition, expControl, tagsStruct, parpoolConn)
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

%% try running simulation
try
       
    
    %% Generate the slice coordinates and transformations
    [spinModel] = domain.generateSimulationDomain( ...
        expControl, acquisition );
    
    %% create experiment progress columns and update to 1%
    exec(expControl.connLocalDB, ...
    ['INSERT INTO edt_tool_local.experiments_progress (exper_id,progress)',...
        ' VALUES (',num2str(expControl.experimentID),',1)']);
    
    %% update progress bar to 1%
    expControl.progress = expControl.progress +1;
    eduTool.frontend.progressUpdate(expControl);
    
    %% NOTE: HERE WE CAN ADD JOB GENERATOR (GPU distribution)
    % Suggestion:   each slice is an independent job
    %               for 2D multi-slice, each slice is simulated and recon
    %               for 3D, still use slices, just sum timeSignals after
    
    %% NOTE: AFTER JOB GENERATOR, PARALLEL EXEC -- JOBS SENT TO QUEUE
    % All below here until recon can be done in parallel
    % For 2D multi-slice, recon can also be done in parallel
    
    %% Get K-space info
    [acquisition] = sequence.tools.generateEncodingInfo(...
        acquisition, expControl );
    
    if expControl.useOldSequence
        %% OLD sequences
        [pulseSequence] = sequence.oldPulseSequenceInterface( ...
            acquisition, mrSystem, expControl );
        % correct number of FE
        acquisition.kSpaceInfo.numFE = 2*acquisition.kSpaceInfo.numFE;
        acquisition.kSpaceInfo.xSize = 2*acquisition.kSpaceInfo.xSize;
        acquisition.kSpaceInfo.feIncidence = 1:acquisition.kSpaceInfo.numFE;
    else
        %% USE NEW CODE
        %% Generate the pulse sequence from the acquisition data
        [pulseSequence] = sequence.generatePulseSequence( ...
            acquisition, mrSystem, expControl );
    end
    
    %% if sequence is PG (diffusion), rotate the diffusion encoding grads
    if isfield(pulseSequence, 'gdwSignal') && ~isempty(pulseSequence.gdwSignal)
        expControl.simulation.simulationEngine = 'diffusion';
        plane = spinModel.slice3D.plane;
        rotMat = (plane.rotMatX*plane.rotMatY*plane.rotMatZ).';
        pulseSequence.gdwSignal = pulseSequence.gdwSignal*rotMat;
    end
    
    %% generate motion sequence
    [motionModel] = motion.generateMotionSequence(pulseSequence, expControl);
    
    %% SAR evaluation, including TIRL
    [sarReport] = coils.evaluateSAR( pulseSequence, ...
        coilSystem, acquisition.data, expControl );
    
    %% update progress bar by 2%
    expControl.progress = expControl.progress +2;
    eduTool.frontend.progressUpdate(expControl);
    
    %% check status
    eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)

    
    %% loop on the number of simulations
    for simNum = 1:expControl.simulation.numberOfSim
        
        %% report start of simulation
        if expControl.debug.debugMode
            tSimulation = tic();
            fprintf(fid, ...
                '\n\n%s : starting simulation %d/%d',...
                functionName, simNum, expControl.simulation.numberOfSim);
            fprintf(fid, '\n');
        end
        
        %% simulate either the 3D model
        if spinModel.is3D
            %% get the 3D domain plane
            plane = spinModel.slice3D.plane;
        else
            %% otherwise, get the plane of the slice number
            plane = spinModel.slice{simNum}.plane;
        end
                    
        %% generate the model by interpolation
        [model] = domain.generateSimulationModel( ...
            plane, anatomicalModel, coilSystem, expControl );

        %% update progress bar by 3%
        expControl.progress = expControl.progress +...
            3/expControl.simulation.numberOfSim;
        eduTool.frontend.progressUpdate(expControl);
        
        %% run the simulation 
        timeSolution = [];
        [timeSolution] = simulator.runSimulation( ...
             model, pulseSequence, motionModel, expControl, timeSolution );
         
        %% update progress bar (simulation assumed 90%)
        expControl.progress = expControl.progress +...
            90/expControl.simulation.numberOfSim;
        eduTool.frontend.progressUpdate(expControl);

         %% check status
         eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)
        
        %% TODO: signal processing -- noise addition, etc
        
        %% perform the reconstruction
        %   assembles K-space and generates the image
        [kSpace, iSpace] = reconstructor.runReconstruction(...
            timeSolution, acquisition, plane, coilSystem, expControl );

        %% copy data into corresponding structures
        if spinModel.is3D
            % if 0, it means that is 3D
            spinModel.slice3D.model = model;
            reconData.slice3D.timeSolution  = timeSolution;
            reconData.slice3D.kSpace        = kSpace;
            reconData.slice3D.iSpace        = iSpace;
        else
            % otherwise, the slice number
            spinModel.slice{simNum}.model = model;
            reconData.slice{simNum}.timeSolution  = timeSolution;
            reconData.slice{simNum}.kSpace        = kSpace;
            reconData.slice{simNum}.iSpace        = iSpace;
        end
        spinModel.totalIsochromats = spinModel.totalIsochromats ...
            + model.numIsochromats;
        
        figure();
        numZ = size(iSpace,3);
        numCI = size(iSpace,5);
        count = 0;
        for ci=1:numCI
            for zz=1:numZ 
                % original image
                count = count+1;
                subplot(numCI,numZ,count);
                iMap = abs(iSpace(:,:,zz,1,ci).');
                imagesc(iMap);
                axis image; colormap gray; colorbar;
                title(sprintf('slice %d contrast %d', zz, ci));
            end
        end
        

        %% report end of experiment
        if expControl.debug.debugMode
            fprintf(fid, ...
                '\n%s : done simulation and reconstruction %d/%d',...
                functionName, simNum, expControl.simulation.numberOfSim);
            fprintf(fid, '\n  Elapsed Time   %.3fs', toc(tSimulation));
            fprintf(fid, '\n\n');
        end
        
    end
    
%     %% process the image and send data back to FE
%     [imageData] = dicom.processReconstruction( ...
%         reconData, acquisition );

%% update progress bar to full
expControl.progress = 100;
eduTool.frontend.progressUpdate(expControl);

catch ME
    
    %% catch the error
    ME.identifier;
    ME.message;
           
    if strcmp(ME.message,'CANCELLED-BY-USER')
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
            expControl.debug.timeStamp, ...
            expControl.userID,...
            expControl.experimentID,...
            ME.stack(1).name,ME.stack(1).line,ME.message );
        fprintf(1, '%s\n', errorMessage);
        
    else % no developer
        
        %% save data for debugging
        errorDataFile = sprintf('%serrorBackend-%s-usr%s-exp%s.mat',...
            expControl.debug.errorFolder, expControl.debug.timeStamp, ...
            expControl.debug.userID, expControl.debug.experimentID);
        save(errorDataFile, 'expControl', 'acquisition', 'ME', '-v7.3');
        
        %% send error in connection to DB for backend
        errorMessage = sprintf(['%s - User (%d) - Exper (%d) - ',...
            'Error in function %s() at line %d.',...
            '\n\nError Message:\n%s',...
            '\n\nData for error replication saved in %s'], ...
            expControl.debug.timeStamp, ...
            expControl.debug.userID,...
            expControl.debug.experimentID,...
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
    
    %% plot sequence
    figure();
    plot(pulseSequence.time,pulseSequence.rfmSignal*1e3);
    xlabel('time (s)');
    hold on
    plot(pulseSequence.time,pulseSequence.gxSignal);
    plot(pulseSequence.time,pulseSequence.gySignal);
    plot(pulseSequence.time,pulseSequence.gzSignal);
    plot(pulseSequence.time(pulseSequence.rxSignal>0), ...
        pulseSequence.gxSignal(pulseSequence.rxSignal>0), 'o');
    plot(pulseSequence.time(pulseSequence.rxLimits(:,1)), ...
        pulseSequence.gxSignal(pulseSequence.rxLimits(:,1)), '>');
    plot(pulseSequence.time(pulseSequence.rxLimits(:,2)), ...
        pulseSequence.gxSignal(pulseSequence.rxLimits(:,2)), '<');
    plot(pulseSequence.time(pulseSequence.partLimits(:,1)),zeros(pulseSequence.numParts,1), '^');
    plot(pulseSequence.time(pulseSequence.partLimits(:,2)),zeros(pulseSequence.numParts,1), 'v');
    if nnz(pulseSequence.swcSignal) > 0
        plot(pulseSequence.time(pulseSequence.swcSignal>0), ...
            zeros(nnz(pulseSequence.swcSignal),1), 's', 'LineWidth', 2, 'MarkerSize', 10);
        legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end', 'SWC');
    else
        legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end');
    end
    
end
