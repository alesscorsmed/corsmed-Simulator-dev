clear all; close all; clc;

corsmedian1 = load('/efs-mount-point/edt_tool/FILES/anatomical_models/20210108_SegarsModel_Temp.mat');

basePath = '/efs-mount-point/S20/INPUTS/anatomical/edutool/CORSMEDIAN/';
fileName = 'corsmedian_20210331_multiPhase.mat';
anatomicalModel = load(sprintf('%s%s',basePath,fileName));

anatomicalModel.tissueValues(:,:) = corsmedian1.model_tissues(:,:);

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


name = 'Varys';
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
