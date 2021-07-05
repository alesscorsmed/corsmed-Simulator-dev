function [spinModel, pulseSequence, motionModel, sarReport, ...
    imageData, reconData, expControl] = slicerV3(...
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
functionName = 'eduTool.run.slicerV3';
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
eduTool.frontend.updateScannerStatus(expControl.connLocalDB, ...
    'Slice and Dice...');

%% set progress bar
expControl.progress = 1;
eduTool.frontend.progressUpdate(expControl);

%% GENERAL: DICOM and DOMAIN DISCRETIZATION

% general dicom info
[dicomInfo] = image.dicom.generateDicomInfo( ...
    acquisition, mrSystem, expControl);

% simulation domain slicing and interpolation
[fovDomain, fovModel] = domain.generateSlices( ...
    acquisition, anatomicalModel, coilSystem, expControl );

%% ACQUISITION RELATED INFO: ENCODING, SEQUENCE, MOTION, SAR

% add an extra field locally to acquisition.data 
% with slice z coordinates after rotation (from fovDomain)
acquisition.data.zSliceCoord = fovDomain.slab.zSliceCoord;

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
        acquisition, encodingData, mrSystem, expControl );
end

%% CORRECTIONS: rotation of the sequence
% if sequence is PG (diffusion), rotate the diffusion encoding grads
if isfield(pulseSequence, 'gdwSignal') && ~isempty(pulseSequence.gdwSignal)
    expControl.simulation.simulationEngine = 'diffusion';
    plane = fovDomain.slab.plane;
    rotMat = (plane.rotMatX*plane.rotMatY*plane.rotMatZ).';
    pulseSequence.gdwSignal = pulseSequence.gdwSignal*rotMat;
end


% generate motion sequence
[motionModel] = motion.generateMotionSequence( pulseSequence, expControl );

% SAR evaluation, including TIRL
[sarReport] = coils.evaluateSAR( pulseSequence, ...
    coilSystem, acquisition.data, expControl );

%% set progress bar
expControl.progress = 2;
eduTool.frontend.progressUpdate(expControl);

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
numGPUs     = expControl.simulation.numGPU;
numSlices   = fovModel.numSlices;
numParts    = lcm(numSlices,numGPUs)/numSlices;
spinModel.numParts  = numParts;
spinModel.numSlices = numSlices;

% loop and assign
jobNum = 0;
for sliceNum = 1:numSlices
    % get interpolated model from the slice
    model = fovModel.slice{sliceNum}.model;
    
    % generate the limits for each partition
    numIso = model.numIsochromats;
    if (numParts > 1) && (numIso > 1024*numParts)
        partLimits = 1:ceil(numIso/numParts):numIso;
        partLimits = [partLimits; [partLimits(2:numParts)+1,numIso] ].';
        localParts = numParts;
    else
        if numIso > 2^18
            numParts = ceil(numIso/2^18);
            partLimits = 1:ceil(numIso/numParts):numIso;
            partLimits = [partLimits; [partLimits(2:numParts)+1,numIso] ].';
            localParts = numParts;
        else
            partLimits = [1, numIso];
            localParts = 1;
        end
    end
    
    for partNum = 1:localParts
        % assign number and corresponding slice and part
        jobNum = jobNum + 1;
        spinModel.data{jobNum}.sliceNum = sliceNum;
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
        partModel.pd = reshape(model.pd(indexes),numIso,1);
        % tissues
        partModel.tissueValues = model.tissueValues;
        partModel.tissueDiff   = model.tissueDiff;
        partModel.tissueType   = reshape(model.tissueType(indexes),numIso,1);
        % coil maps
        partModel.rxCoilMapsX = reshape(model.rxCoilMapsX(indexes,:),numIso,[]);
        partModel.rxCoilMapsY = reshape(model.rxCoilMapsY(indexes,:),numIso,[]);
        partModel.numRxCoils  = model.numRxCoils;
        % assign
        spinModel.data{jobNum}.model = partModel;
        spinModel.data{jobNum}.GPUindex = eduTool.multiGPU.chooseGPU(jobNum);
    end
end
% update the actual number of jobs
spinModel.numJobs = jobNum;
expControl.simulation.numberOfSim = jobNum;

%% update progress bar to 5%
expControl.progress = 5;
eduTool.frontend.progressUpdate(expControl);

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Number of Jobs    %d', spinModel.numJobs);
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
