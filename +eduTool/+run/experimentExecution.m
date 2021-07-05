function [stats, testControl] = experimentExecution(...
    anatomicalModel, coilSystem, acquisition, expControl, sessionData, testControl)
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
if (nargin < 4)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end
if (nargin < 6) || isempty(testControl)
    testControl.testGen = 0;
    testControl.testRun = 0;
end

% %% return for MOLLI
% if strcmpi(acquisition.data.pulseSeqFamilyName, 'molli')
%     msg = sprintf( ['The selected sequence type (%s) is not available in V2. ',...
%         'We will upgrade V2 with the required sequence in the coming release. '...
%         'In the meantime please switch to V1: ' ...
%         ' Simulation -> Advaced -> EduTool version -> V1 ', ...
%         ' and  Simulation -> Advaced -> Pulse seq. generator -> V1'],...
%         acquisition.data.pulseSeqFamilyName );
%     eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
% end


%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    fprintf(fid, '\n%s : start', functionName);
end
% time it
startTStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
tTotal = tic();

%% SLICER
tSlicer = tic();
[spinModel, pulseSequence, sarReport, imageData, reconData, ...
    expControl] = eduTool.run.slicer( anatomicalModel, coilSystem, ...
    expControl.mrSystem, acquisition, expControl );
tSlicer = toc(tSlicer);

%% verify testing
if testControl.testRun
    %% we are running testing
    [testControl] = eduTool.test.regressionVerifyStruct( spinModel, testControl, 'spinModel');
    [testControl] = eduTool.test.regressionVerifyStruct( pulseSequence, testControl, 'pulseSequence');
    [testControl] = eduTool.test.regressionVerifyStruct( sarReport, testControl, 'sarReport');
end
if testControl.testGen
    %% we are generating test data
    [testControl] = eduTool.test.regressionSaveStruct( spinModel, testControl, 'spinModel');
    [testControl] = eduTool.test.regressionSaveStruct( pulseSequence, testControl, 'pulseSequence');
    [testControl] = eduTool.test.regressionSaveStruct( sarReport, testControl, 'sarReport');    
end

%% check status
eduTool.frontend.checkExperimentStatus(expControl,expControl.experimentID);

%% SIMULATOR
simStartTStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
tEngine = tic();
if isfield(expControl.simulation, 'parPool') ...
        && (expControl.simulation.parPool.NumWorkers > 1)
    [simSignal] = eduTool.run.engineParallel( ...
        spinModel, pulseSequence, expControl);
else
    [simSignal] = eduTool.run.engine( ...
        spinModel, pulseSequence, expControl );
end
tEngine = toc(tEngine);
simFinalTStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');

%% verify testing
if testControl.testRun
    %% we are running testing
    [testControl] = eduTool.test.regressionVerifyDataStruct( simSignal, testControl, 'simSignal');
end
if testControl.testGen
    %% we are generating test data
    [testControl] = eduTool.test.regressionSaveStruct( simSignal, testControl, 'simSignal');  
end

%% check status
eduTool.frontend.checkExperimentStatus(expControl,expControl.experimentID)

%% Make the simulator to crash for test purposes - REMOVE IT @@@ CX
fprintf(1,'\nDEBUGGING:\n')
fprintf(1,'Approach: %s\n',expControl.approach)
fprintf(1,'Pulseq family: %s\n',acquisition.data.pulseSeqFamilyName)
fprintf(1,'Slice thickn.: %s\n',num2str(acquisition.data.sliceThickness))

if strcmp(expControl.approach,'jsonstandalone')
    if strcmpi(acquisition.data.pulseSeqFamilyName, 'ss-fse')
        if strcmp(num2str(acquisition.data.sliceThickness),num2str(0.007))
            ME = MException('eduTool:errorOnPurpose',...
                '%s : error for test purposes',functionName);
            throw(ME);
        elseif strcmp(num2str(acquisition.data.sliceThickness),num2str(0.008))
            msg = sprintf('Corsmed error for test purposes (CANCELLED-ERROR)!');
            eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
        elseif strcmp(num2str(acquisition.data.sliceThickness),num2str(0.009))
            msg = sprintf('Corsmed error for test purposes (CONFIRM)!');
            eduTool.frontend.askUserToConfirmExperiment(expControl,sessionData,msg);            
        end
    end
end

%% RECONSTRUCTION
tRecon = tic();
[reconData] = eduTool.run.recon(...
    simSignal, reconData, expControl );
tRecon = toc(tRecon);

%% verify testing
if testControl.testRun
    %% we are running testing
    [testControl] = eduTool.test.regressionVerifyDataStruct( reconData, testControl, 'reconData');
end
if testControl.testGen
    %% we are generating test data
    [testControl] = eduTool.test.regressionSaveStruct( reconData, testControl, 'reconData');  
end

%% check status
eduTool.frontend.checkExperimentStatus(expControl,expControl.experimentID)

%% DICOM GENERATION
if isfield(expControl,'mode') && strcmp(expControl.mode,'test')
    tDicom = 0;
else
    tDicom = tic();
    [imageData,~,jsonStructure] = eduTool.run.dicom(...
        reconData, imageData, expControl, sarReport);
    tDicom = toc(tDicom);
end

%% verify testing
if testControl.testRun
    %% mask air to remove phase and other issues
    backgroundTissue = anatomicalModel.backgroundTissue;
    for ss = 1:imageData.numSlices
        % small data image intesity
        for cc = 1:imageData.numContrasts
            idxZero = imageData.slice{ss}.contrast{cc}.mask == backgroundTissue;
            imageData.slice{ss}.contrast{cc}.image(idxZero) = 0;
            imageData.slice{ss}.contrast{cc}.normalizedImage(idxZero) = 0;
            if ( strcmpi(imageData.processContrast,'waterfatoop') && imageData.b0map == 1) ...
                    || strcmpi(imageData.processContrast,'swi')
                % Processed  and FFT resulting images are nuts for phase
                imageData.slice{ss}.contrast{cc}.normalizedImage(:) = 0;
                imageData.slice{ss}.contrast{cc}.kspaceOutputComplex(:) = 0;
                imageData.slice{ss}.contrast{cc}.kspaceOutputImage(:) = 0;
            end
        end
    end
    %% we are running testing
    [testControl] = eduTool.test.regressionVerifyDataStruct( imageData, testControl, 'imageData');
end
if testControl.testGen
    %% handle angle maps
    backgroundTissue = anatomicalModel.backgroundTissue;
    for ss = 1:imageData.numSlices
        % small data image intesity
        for cc = 1:imageData.numContrasts
            idxZero = imageData.slice{ss}.contrast{cc}.mask == backgroundTissue;
            imageData.slice{ss}.contrast{cc}.image(idxZero) = 0;
            imageData.slice{ss}.contrast{cc}.normalizedImage(idxZero) = 0;
            if ( strcmpi(imageData.processContrast,'waterfatoop') && imageData.b0map == 1) ...
                    || strcmpi(imageData.processContrast,'swi')
                % Processed  and FFT resulting images are nuts for phase
                imageData.slice{ss}.contrast{cc}.normalizedImage(:) = 0;
                imageData.slice{ss}.contrast{cc}.kspaceOutputComplex(:) = 0;
                imageData.slice{ss}.contrast{cc}.kspaceOutputImage(:) = 0;
            end
        end
    end
    %% we are generating test data
    [testControl] = eduTool.test.regressionSaveStruct( imageData, testControl, 'imageData');  
end

% total time
tTotal = toc(tTotal);
finalTStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');

% include exectime in jsonStructure for use in DB
if strcmpi(expControl.application,'edutool') && ...
    isfield(expControl,'approach') && ...
    strcmp(expControl.approach,'jsonstandalone')

    %execution time includes in UIindicators struct
    jsonStructure.UIindicators.timeExec = tTotal;

    expControl.redis.R = tools.redis.redisSetJsonWrapper(expControl.redis.R,...
    expControl.redis.keys.experimentResultsRedisKey,jsonStructure);
else
    eduTool.frontend.updateExecTime(expControl,tTotal)
end

%% update progress bar to full
expControl.progress = 100;
eduTool.frontend.updateExperimentProgress(expControl);

%% assign stats
stats.userID         = expControl.userID;
stats.courseID       = expControl.courseID;
stats.experimentID   = expControl.experimentID;
stats.expName        = expControl.name;
stats.startTStamp    = startTStamp;
stats.finalTStamp    = finalTStamp;
stats.simStartTStamp = simStartTStamp;
stats.simFinalTStamp = simFinalTStamp;

stats.txCoil        = coilSystem.coilModel{coilSystem.indexTx}.data.name;
stats.rxCoil        = coilSystem.coilModel{coilSystem.indexRx}.data.name;
stats.numRxCoils    = coilSystem.coilModel{coilSystem.indexRx}.data.numCoils;
stats.anatomy       = imageData.name;
stats.bodyPart      = imageData.bodyPartName;

stats.seqName       = pulseSequence.name;
stats.seqFamily     = pulseSequence.familyName;
stats.seqTIRL       = pulseSequence.totalTime;
stats.numSteps      = pulseSequence.numSteps;
stats.numRxs        = pulseSequence.numRxs;
stats.numVoxels     = spinModel.totalIso;
stats.numSlices     = spinModel.numSlices;
stats.numJobs       = spinModel.numJobs;

stats.matrixSizeX    = acquisition.data.numFE;
stats.matrixSizeY    = acquisition.data.numPE;

stats.voxelSizeX     = acquisition.data.fovFE*1000/acquisition.data.matrixX; 
stats.voxelSizeY     = acquisition.data.fovPE*1000/acquisition.data.matrixY; 
stats.voxelSizeZ     = acquisition.data.sliceThickness*1000;


stats.typeGPU       = expControl.simulation.gpuName;
stats.numCPUs       = expControl.simulation.numCPUs;
stats.numGPUs       = expControl.simulation.numGPUs;
stats.numThreads    = expControl.simulation.threads;
stats.kernel        = expControl.simulation.kernelPtx;

stats.timeEstRun    = expControl.estimatedRunTime;
stats.timeSlicer    = tSlicer;
stats.timeEngine    = tEngine;
stats.timeRecon     = tRecon;
stats.timeDicom     = tDicom;
stats.timeTotal     = tTotal;

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for experiment %d, elapsed time %.3fs',...
        functionName, expControl.experimentID, tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
