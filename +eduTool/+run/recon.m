function [reconData] = recon( simSignal, reconData, expControl )
%
% EDUTOOL.RUN.RECON
%
%	Runs the reconstruction.
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.recon';
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

%% if there is no DB connection
if ~isfield(expControl,'connLocalDB')
    expControl.connLocalDB = [];
end

%% RECONSTRUCTION INSTANCE
eduTool.frontend.updateScannerStatus(expControl,'Signal to Visual...');

%% prepare common data
numSlices = reconData.numSlices;
numFrames = reconData.numFrames;
numCoils  = simSignal.numCoils;
numReads  = simSignal.numReads;
% allocate, zero signal data, keep time
rxTime = simSignal.timeSolution{1}.time;
    
%% allocate space in each frame and slice
for frameIdx = 1:numFrames
    for sliceIdx = 1:numSlices
        reconData.slice{sliceIdx}.frame{frameIdx}.raw.numCoils  = numCoils;
        reconData.slice{sliceIdx}.frame{frameIdx}.raw.numReads  = numReads;
        reconData.slice{sliceIdx}.frame{frameIdx}.raw.time      = rxTime;
        reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sx        = zeros(numReads,numCoils);
        reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sy        = zeros(numReads,numCoils);
        reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sz        = zeros(numReads,numCoils);
    end
end

%% combine all simulated data in correct k-spaces
for simNum = 1:simSignal.numJobs
    % get the slice
    sliceIdx = simSignal.timeSolution{simNum}.sliceNum;
    %phaseIdx = simSignal.timeSolution{simNum}.phaseNum;
    %contrIdx = simSignal.timeSolution{simNum}.contrNum;
    frameIdx = simSignal.timeSolution{simNum}.frameNum;
    % add signals to the slice
    reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sx = ...
        + reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sx ...
        + simSignal.timeSolution{simNum}.Sx;
    reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sy = ...
        + reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sy ...
        + simSignal.timeSolution{simNum}.Sy;
    reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sz = ...
        + reconData.slice{sliceIdx}.frame{frameIdx}.raw.Sz ...
        + simSignal.timeSolution{simNum}.Sz;
end
    
%% loop on the slices and recon
for frameIdx = 1:numFrames
    for sliceIdx = 1:numSlices
        
        %% extract raw struct with signals
        raw = reconData.slice{sliceIdx}.frame{frameIdx}.raw;
        
        %% noise generation: update raw struct with noise signal field
        [raw] = noise.generateNoiseSignal(raw,...
            reconData.noise, expControl);
        
        %% perform the reconstruction
        %   assembles K-space and generates the image
        [kSpace, iSpace] = reconstructor.runReconstruction(...
            raw, reconData.slice{sliceIdx}.sens,...
            reconData.encoding, expControl );
        
        %% copy recon data into corresponding structures
        reconData.slice{sliceIdx}.frame{frameIdx}.raw    = raw;
        reconData.slice{sliceIdx}.frame{frameIdx}.kSpace = kSpace;
        reconData.slice{sliceIdx}.frame{frameIdx}.iSpace = iSpace;
    end
end

%% other info
[numX,numY,numZ,~,numCI] = size(iSpace);
reconData.numX = numX;
reconData.numY = numY;
reconData.numZ = numZ;
reconData.numC = numCI;

%% update progress bar to 99
expControl.progress = 99;
eduTool.frontend.updateExperimentProgress(expControl);

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Number of Slices  %d', reconData.numSlices);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end