function [spinModel, pulseSequence, sarReport, ...
    imageData, reconData, expControl] = slicer(...
    anatomicalModel, coilSystem, mrSystem, acquisition, expControl)
%
% EDUTOOL.RUN.SLICER
%
%	Slice and Dice.
%
% INPUT
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.slicer';
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

%% if there is no DB connection
if ~isfield(expControl,'connLocalDB')
    expControl.connLocalDB = [];
end

%% SLICER INSTANCE
eduTool.frontend.updateScannerStatus(expControl,'Slice and Dice...');

%% set progress bar
expControl.progress = 1;
eduTool.frontend.updateExperimentProgress(expControl);

%% GENERAL: DICOM and DOMAIN DISCRETIZATION

% general dicom info
[dicomInfo] = image.dicom.generateDicomInfo( ...
    acquisition, mrSystem, expControl);

% simulation domain slicing and plane generation
[fovDomain] = domain.generateSimulationDomain( acquisition, expControl );

%% ACQUISITION RELATED INFO: ENCODING, SEQUENCE, MOTION, SAR

% add an extra field locally to acquisition.data 
% with slice z coordinates after rotation (from fovDomain)
acquisition.data.zSliceCoord = fovDomain.slab.zSliceCoord;

% if (contains(lower(acquisition.data.pulseSeqFamilyName), 'gre-oop') || ...
%         contains(lower(acquisition.data.pulseSeqFamilyName), 'swi-2d-gre'))
%     acquisition.data.samplingFactorFE = 1;
% end

% Generate the encoding information (encoding planner)
[encodingData] = encoder.generateEncodingPlan( acquisition, expControl );

% Generate the sequence
if expControl.useOldSequence
    % OLD sequences
    [pulseSequence, encodingData] = sequence.oldPulseSequenceInterface( ...
        acquisition, encodingData, mrSystem, expControl );
else
    % USE NEW CODE
    % Generate the pulse sequence from the acquisition data
    [pulseSequence] = sequence.generatePulseSequence( ...
        acquisition, encodingData, mrSystem, expControl, anatomicalModel);
end

%% CORRECTIONS: rotation of the sequence
% if sequence is PG (diffusion), rotate the diffusion encoding grads
if isfield(pulseSequence, 'gdwSignal') && ~isempty(pulseSequence.gdwSignal)
    expControl.simulation.simulationEngine = 'diffusion';
    plane = fovDomain.slab.plane;
    rotMat = (plane.rotMatX*plane.rotMatY*plane.rotMatZ).';
    pulseSequence.gdwSignal = pulseSequence.gdwSignal*rotMat;
end

% Apply perfusion
if contains(lower(acquisition.data.pulseSeqFamilyName), 'perfusion')
    expControl.model.perfusion.apply = 1;
end

% SAR evaluation, including TIRL, and front end update
[sarReport] = coils.evaluateSAR( pulseSequence, ...
    coilSystem, acquisition.data, expControl, anatomicalModel );

%% set progress bar
expControl.progress = 2;
eduTool.frontend.updateExperimentProgress(expControl);

%% GENERATE THE INTERPOLATED SLICES
% model interpolation
[fovModel] = simModel.generateSimulationModel( ...
    fovDomain, anatomicalModel, coilSystem, ...
    expControl.model, expControl.debug );

% compute total number of frames: phases x contrasts
numSlices = fovModel.numSlices;
if strcmp(acquisition.data.pulseSeqFamilyName,'cine-bssfp')
    % set frames to be passed to simulation
    [phasesIdx, numPhases]          = sequence.tools.CINEcalculatePhases(...
        acquisition,anatomicalModel);
    numContrasts                    = 1;
    numFrames                       = numPhases*numContrasts;
    encodingData.plan.numPhases     = numPhases;
    encodingData.plan.phaseIndex    = rem(phasesIdx,fovModel.numPhases+1); % array of phases to simulate
    encodingData.plan.numContrasts  = numContrasts;
    encodingData.plan.numFrames     = numFrames;
else
    numPhases                       = 1;
    if expControl.model.perfusion.apply
        numContrasts = expControl.model.perfusion.contrasts;
    else
        numContrasts                = 1;
    end
    numFrames                       = numPhases*numContrasts;
    encodingData.plan.numPhases     = numPhases;
    encodingData.plan.phaseIndex    = expControl.model.cardiacPhase; % use the default cardiac phase
    encodingData.plan.numContrasts  = numContrasts;
    encodingData.plan.numFrames     = numFrames; 
end


%% RECON RELATED INFO: IMAGE and RECON

% image slices and data
[imageData] = image.generateSlices( ...
    fovDomain, pulseSequence, encodingData.plan, expControl );

% upgrade imageData with dicomInfo
imageData.dicomInfo = dicomInfo;

% reconstruction slices and data
[reconData] = reconstructor.generateSlices( ...
    fovDomain, anatomicalModel, coilSystem, ...
    encodingData, acquisition.noise, expControl );


%% set progress bar
expControl.progress = 4;
eduTool.frontend.progressUpdate(expControl);

%% SIMULATION RELATED INFO: GENERATE SIMULATION PARTS

% number of parts to get the least common multiple of slices and GPUs
numGPUs     = expControl.simulation.numGPUs;
numCPUs     = expControl.simulation.numCPUs;
numPROC     = max(numGPUs,numCPUs);
numParts    = lcm(numSlices,numPROC)/numSlices;
spinModel.numParts  = numParts;
spinModel.numSlices = numSlices;

% extract info to dermine partition size:
estPerf     = 50e-12; % time / voxel / timeStep
numSteps    = pulseSequence.numSteps;
avgPartTime = 5; % partition so we estimated each job to take 5 sec
maxPartSize = floor(log2(avgPartTime/(numSteps*estPerf)));
maxPartSize = max(2^maxPartSize,2^16);

% loop and assign
jobNum = 0;
totalIso = 0;
for sliceNum = 1:numSlices
    
    % get interpolated model from the slice
    model = fovModel.slice{sliceNum}.model;
    
    % generate the limits for each partition
    numIso = model.numIsochromats;
    localParts = numParts;
    % decide if model is too large, and we require further partition
    if numIso/localParts > maxPartSize
        localParts = numPROC*ceil(numIso/(numPROC*maxPartSize));
        partLimits = 1:ceil(numIso/localParts):numIso;
        partLimits = [partLimits; [partLimits(2:localParts)+1,numIso] ].';
    else
        % other wise, check if partition too small
        if (localParts > 1) && (numIso > 1024*localParts)
            partLimits = 1:ceil(numIso/localParts):numIso;
            partLimits = [partLimits; [partLimits(2:localParts)+1,numIso] ].';
        else
            % if too small, single part
            partLimits = [1, numIso];
            localParts = 1;
        end
    end
    
    % generate the individual models
    for contrNum = 1:numContrasts
        for phaseNum = 1:numPhases
            for partNum = 1:localParts
                % assign number and corresponding slice and part
                jobNum = jobNum + 1;
                spinModel.data{jobNum}.sliceNum = sliceNum;
                spinModel.data{jobNum}.phaseNum = phaseNum;
                spinModel.data{jobNum}.contrNum = contrNum;
                spinModel.data{jobNum}.frameNum = (contrNum-1)*numPhases + phaseNum;
                spinModel.data{jobNum}.numParts = localParts;
                spinModel.data{jobNum}.partNum  = partNum;
                % assig model data
                indexes = partLimits(partNum,1):partLimits(partNum,2);
                numIso = length(indexes);
                % size and resolution
                partModel.numIsochromats = numIso;
                partModel.resolution     = model.resolution;
                % constants
                partModel.mu = model.mu;
                partModel.b0 = model.b0;
                % coordinates
                partModel.x  = reshape(model.x(indexes),numIso,1);
                partModel.y  = reshape(model.y(indexes),numIso,1);
                partModel.z  = reshape(model.z(indexes),numIso,1);
                % voxel properties
                partModel.bi = reshape(model.bi(indexes),numIso,1);
                if isempty(model.pd)
                    partModel.pd = [];
                else
                    partModel.pd = reshape(model.pd(indexes),numIso,1);
                end
                % tissues: use contrast and frame
                if expControl.model.perfusion.apply
                    % apply perfusion
                    partModel.tissueValues = ...
                        simModel.perfusion.perfusionUpdateTissueProp( ...
                        model.tissueValues, expControl.model.perfusion, contrNum); %% TODO: assign given contrast
                else
                    % unchanged tissue  values
                     partModel.tissueValues = model.tissueValues;
                end
                % apply the corresponding phase
                phaseIndex = encodingData.plan.phaseIndex(phaseNum); % select the correspondig phase number from the array
                partModel.tissueType   = reshape(model.tissueType(indexes,phaseIndex),numIso,1);
                partModel.tissueDiff   = model.tissueDiff;
                % coil maps
                partModel.rxCoilMapsX = reshape(model.rxCoilMapsX(indexes,:),numIso,[]);
                partModel.rxCoilMapsY = reshape(model.rxCoilMapsY(indexes,:),numIso,[]);
                partModel.numRxCoils  = model.numRxCoils;
                % assign
                spinModel.data{jobNum}.model = partModel;
                spinModel.data{jobNum}.GPUindex = eduTool.multiGPU.chooseGPU(jobNum);
                totalIso = totalIso + numIso; % acumulate total isochromats
            end
        end
    end
end
% update the actual number of jobs
spinModel.numJobs   = jobNum;
spinModel.totalIso  = totalIso;
expControl.simulation.numberOfSim = jobNum;

%% update progress bar to 5%
expControl.progress = 5;
eduTool.frontend.progressUpdate(expControl);

%% estimate run time based on steps and isochromats
estRunTime = (pulseSequence.numSteps*totalIso)*estPerf;
estRunTime = estRunTime + jobNum*1; % add set up time for sim (GPU trans, so on)
expControl.estimatedRunTime = estRunTime;

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time        %.3fs', tTotal);
    fprintf(fid, '\n  Number of Jobs      %d', spinModel.numJobs);
    fprintf(fid, '\n  Estimated run time  %.3fs', estRunTime);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
    %     %% plot sequence
    %     figure();
    %     plot(pulseSequence.time,pulseSequence.rfmSignal*1e3);
    %     xlabel('time (s)');
    %     hold on
    %     plot(pulseSequence.time,pulseSequence.gxSignal);
    %     plot(pulseSequence.time,pulseSequence.gySignal);
    %     plot(pulseSequence.time,pulseSequence.gzSignal);
    %     plot(pulseSequence.time(pulseSequence.rxSignal>0), ...
    %         pulseSequence.gxSignal(pulseSequence.rxSignal>0), 'o');
    %     plot(pulseSequence.time(pulseSequence.rxLimits(:,1)), ...
    %         pulseSequence.gxSignal(pulseSequence.rxLimits(:,1)), '>');
    %     plot(pulseSequence.time(pulseSequence.rxLimits(:,2)), ...
    %         pulseSequence.gxSignal(pulseSequence.rxLimits(:,2)), '<');
    %     plot(pulseSequence.time(pulseSequence.partLimits(:,1)),zeros(pulseSequence.numParts,1), '^');
    %     plot(pulseSequence.time(pulseSequence.partLimits(:,2)),zeros(pulseSequence.numParts,1), 'v');
    %     if nnz(pulseSequence.swcSignal) > 0
    %         plot(pulseSequence.time(pulseSequence.swcSignal>0), ...
    %             zeros(nnz(pulseSequence.swcSignal),1), 's', 'LineWidth', 2, 'MarkerSize', 10);
    %         legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end', 'SWC');
    %     else
    %         legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end');
    %     end
    
end
