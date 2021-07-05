function [batchName] = generateExperimentName(acqData,expControl)
%
% DOMAIN.GENERATEEXPERIMENTNAME
%
%
% INPUT
%   acquisition        structure with acquisition data
%   expControl         experiment control data, with commLocalDB
%
% OUTPUT
%   batch              name
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.batch.generateExperimentName';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% sequence name
if acqData.is3D
    seqName = sprintf('%s_3D_%dSlices', ...
        acqData.pulseSeqFamilyName,acqData.numSlices);
else
    seqName = sprintf('%s_2D_%dSlices', ...
        acqData.pulseSeqFamilyName,acqData.numSlices);
end
if ~strcmpi(acqData.parallelImaging, 'no')
    seqName = sprintf('%s_sense%d', seqName, acqData.rFactor);
end
switch lower(acqData.partialFourier)
    case lower('readConjugate')
        seqName = sprintf('%s_%s%d', seqName,...
            acqData.partialFourier,round(acqData.fFactor*100));
    case lower('phaseConjugate')
        seqName = sprintf('%s_%s%d', seqName,...
            acqData.partialFourier,round(acqData.fFactor*100));
    case lower('sliceConjugate')
        seqName = sprintf('%s_%s%d', seqName,...
            acqData.partialFourier,round(acqData.fFactor*100));
    otherwise
        seqName = seqName;
end

%% find plane orientation
pointsAll           = acqData.pointsAll;
pointsAllFrontend   = acqData.pointsAllFrontend;  
% extract coordinates data
points          = pointsAll{1};
points          = strsplit(points,';');
pointsFrontEnd  = pointsAllFrontend{1};
pointsFrontEnd  = strsplit(pointsFrontEnd,';');
% generate plane
% assign points
plane.LTop = str2num(points{1,1});  % Top-Left
plane.RTop = str2num(points{1,2});  % Top-Right
plane.LBot = str2num(points{1,3});  % Bottom-Left
plane.RBot = str2num(points{1,4});  % Bottom-right
% get the points from the Front End
plane.LTopFrontEnd  = str2num(pointsFrontEnd{1,1});  % Top-Left
plane.RTopFrontEnd  = str2num(pointsFrontEnd{1,2});  % Top-Right
plane.LBotFrontEnd  = str2num(pointsFrontEnd{1,3});  % Bottom-Left
plane.RBotFrontEnd  = str2num(pointsFrontEnd{1,4});  % Bottom-right
% compute the transformations and final plane coordinates
[plane] = domain.planeHandling.generatePlane(plane);
planeType = strrep(plane.Type, '-','');
if strcmpi(planeType, 'inplane')
    planeType = [planeType,'_',plane.imType];
end

%% motion info
motionSpecs = expControl.motionSpecs;
switch lower(motionSpecs.pattern)
    case 'rotational'
        motionType = 'RotMotion';
    case 'translational'
        switch motionSpecs.transAxis
            case 1 % X
                motionType = 'xTransMotion';
            case 2 % Y
                motionType = 'yTransMotion';
            case 3 % Z
                motionType = 'zTransMotion';
            otherwise % no movement
                motionType = [];
        end
    otherwise
        motionType = [];
end

%% experiment name
if isempty(motionType)
    batchName = sprintf('%s_%s', seqName, planeType);
else
    batchName = sprintf('%s_%s_%s', seqName, planeType, motionType);
end