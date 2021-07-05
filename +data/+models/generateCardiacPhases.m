clear all; close all; clc;
% basePath = '/efs-mount-point/edt_tool/FILES/anatomical_models/CORSMEDIAN/';
basePath = '/efs-mount-point/S20/INPUTS/anatomical/edutool/CORSMEDIAN/';

%% generate the cardiac phases
numFrames   = 20;
bpm         = 60;
cycleTime   = 60/bpm;
frameTime   = linspace(0,cycleTime,numFrames);

% find voxels changing tissue with frame
frameIdx = [];
for ii = 1:numFrames
    fileName = sprintf('corsmedian_20210318_phase%d_v2.mat',ii-1);
    refModel = load(sprintf('%s%s',basePath,fileName));
    for jj = ii+1:numFrames
        fileName = sprintf('corsmedian_20210318_phase%d_v2.mat',jj-1);
        localModel = load(sprintf('%s%s',basePath,fileName));
        [idxRow,idxCol] = find( abs(refModel.tissueType - localModel.tissueType) > 0.0 );
        frameIdx = union(frameIdx, idxRow);
        fprintf(1, '\n Reference frame %d with frame %d processed', ii, jj);
    end
    fprintf(1, '\n');
end
clear refModel; clear localModel;
numVoxelsFrame  = numel(frameIdx);

% create the tissueType of the changing voxels
fileName = 'corsmedian_20210318_phase0_v2.mat';
anatomicalModel = load(sprintf('%s%s',basePath,fileName));
% correct to 6 properties
if anatomicalModel.numProperties < 6
    anatomicalModel.numProperties = 6;
    anatomicalModel.tissueValues(anatomicalModel.numTissues,6) = 0.0;
end
%%TODO
% there are some T1/T2 NANs in the model: set to small T1/T2 values and 0 PD
[idxR,idxC] = find(isnan(anatomicalModel.tissueValues));
anatomicalModel.tissueValues(idxR,1) = 1e-5;
anatomicalModel.tissueValues(idxR,2) = 1e-4;
anatomicalModel.tissueValues(idxR,3) = 0.0;

% assign frame info
anatomicalModel.numFrames   = numFrames;
anatomicalModel.bpm         = bpm;
anatomicalModel.cycleTime   = cycleTime;
anatomicalModel.frameTime   = frameTime;

% tissue types for frame voxels
anatomicalModel.idxFrame        = frameIdx;
anatomicalModel.numVoxelsFrame  = numVoxelsFrame;
anatomicalModel.tissueTypeFrame = zeros(numVoxelsFrame,numFrames);

% fill data for current frame
anatomicalModel.tissueTypeFrame(:,1) = anatomicalModel.tissueType(anatomicalModel.idxFrame,1);
% fill rest of the frames
for ii = 2:numFrames
    fileName = sprintf('corsmedian_20210318_phase%d_v2.mat',ii-1);
    localModel = load(sprintf('%s%s',basePath,fileName));
    anatomicalModel.tissueTypeFrame(:,ii) = localModel.tissueType(anatomicalModel.idxFrame,1);
    fprintf(1, '\n Frame %d data processed', ii);
end
fprintf(1, '\n');

name = 'varys_20210406_multiPhase';
b0Inhomogeneity = anatomicalModel.b0Inhomogeneity;
backgroundTissue = anatomicalModel.backgroundTissue;
dimensions = anatomicalModel.dimensions;
domain = anatomicalModel.domain;
fatTissuesIDs = anatomicalModel.fatTissuesIDs;
isGridded = anatomicalModel.isGridded;
numProperties = anatomicalModel.numProperties;
numTissues = anatomicalModel.numTissues;
pdFluctuationFactor = anatomicalModel.pdFluctuationFactor;
pdInhomogeneity = [];
resolution = anatomicalModel.resolution;
spatial = anatomicalModel.spatial;
tissueType = anatomicalModel.tissueType;
tissueValues = anatomicalModel.tissueValues;

numFrames = anatomicalModel.numFrames;
bpm = anatomicalModel.bpm;
cycleTime = anatomicalModel.cycleTime;
frameTime = anatomicalModel.frameTime;
idxFrame = anatomicalModel.idxFrame;
numVoxelsFrame = anatomicalModel.numVoxelsFrame;
tissueTypeFrame = anatomicalModel.tissueTypeFrame;

fileName = 'varys_20210406_multiPhase.mat';
save(sprintf('%s%s',basePath,fileName),...
'b0Inhomogeneity',...
'backgroundTissue',...
'dimensions',...
'domain',...
'fatTissuesIDs',...
'isGridded',...
'name',...
'numProperties',...
'numTissues',...
'pdFluctuationFactor',...
'pdInhomogeneity',...
'resolution',...
'spatial',...
'tissueType',...
'tissueValues',...
'numFrames',...
'bpm',...
'cycleTime',...
'frameTime',...
'idxFrame',...
'numVoxelsFrame',...
'tissueTypeFrame',...
'-v7.3');


basePath = '/efs-mount-point/S20/INPUTS/anatomical/edutool/CORSMEDIAN/';
fileName = 'corsmedian_20210331_multiPhase.mat';
anatomicalModel = load(sprintf('%s%s',basePath,fileName));