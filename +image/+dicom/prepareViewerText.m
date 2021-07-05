function [textForImage] = prepareViewerText( ...
        imagePlane, imageData, pulseSequence, encodingPlan,...
        expControl, sliceNum, contrastNum )
%
% IMAGE.DICOM.PREPAREVIEWERTEXT
%
%     Creates text for the image viewer
%
% INPUT
%   
%
% OUTPUT
%   
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'image.dicom.prepareViewerText';
if (nargin < 7)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end


%% Top Left: common line
if isfield(pulseSequence,'displaySeqName')
    textForImage.LTop{1} = [num2str(pulseSequence.seqNum),...
        ' - ',pulseSequence.displaySeqName];
else
    textForImage.LTop{1} = [num2str(pulseSequence.seqNum),...
        ' - ',pulseSequence.familyName];
end
if expControl.useOldSequence
    textForImage.LTop{1} = [textForImage.LTop{1},' v1'];
end
textForImage.LTop{2} = [ 'Slice ', ...
    num2str(sliceNum),'/', num2str(imageData.numSlices) ];
lineCount = 3;
if imageData.numContrasts > 1
    textForImage.LTop{lineCount} = [ 'Contrast ', ...
        num2str(contrastNum),'/',num2str(imageData.numContrasts) ];
    lineCount = lineCount + 1;
end
% DW information
if isfield(pulseSequence,'encPG') ...
        && ~isempty(pulseSequence.encPG) ...
        && (pulseSequence.encPG.Beta > 0)
    BETA = pulseSequence.encPG.Beta/1e6; % in s/mm2
    betaScale = floor(log10(BETA));
    textForImage.LTop{lineCount} = sprintf('B (%s): %#1.2fe%d s/mm2',...
    upper(pulseSequence.encPG.Dir), BETA/(10^betaScale), betaScale);
end

%% Top Right
gridStep = expControl.model.gridStep;
textForImage.RTop{1} = imageData.bodyPartName;
textForImage.RTop{2} = ['Grid: ', num2str(gridStep(1)*1000),...
    'x',num2str(gridStep(2)*1000),...
    'x',num2str(gridStep(3)*1000),' mm'];
textForImage.RTop{3} = expControl.model.coilType;

%% Bottom Right
textForImage.RBot{3} = expControl.timeStamp;

%% Bottom Left
% Sequence properties
if contains(lower(pulseSequence.familyName),'epi')
    sequenceText = sprintf('TR/effTE/ESP: %.2f/%.2f/%.2f ms',...
        pulseSequence.TR*1000,...
        pulseSequence.effTE*1000,...
        pulseSequence.ESP*1000);
else
    if pulseSequence.numTL > 1
        sequenceText = sprintf('TR/effTE: %.2f/%.2f ms - ETL: %d',...
            pulseSequence.TR*1000,...
            pulseSequence.effTE*1000,...
            pulseSequence.numTL);
    else
        sequenceText = sprintf('TR/TE: %.2f/%.2f ms',...
            pulseSequence.TR*1000,...
            pulseSequence.TE*1000);
    end
end
% IR
if pulseSequence.TI > 0
    if strcmpi(pulseSequence.familyName,'molli') 
        if contrastNum<=size(pulseSequence.MOLLI.TIs,2)
            irSenseText = sprintf('TI= %.1fms',...
                pulseSequence.MOLLI.TIs(contrastNum)*1e3);
        else
            irSenseText = [];
        end
    else
        irSenseText = sprintf('TI= %.1fms',pulseSequence.TI*1e3);
    end
else
    irSenseText = [];
end
% SENSE
if encodingPlan.rFactorPE > 1
    irSenseText = [irSenseText, sprintf('SENSE: %d', encodingPlan.rFactorPE)];
end
% partial fourier
if ( encodingPlan.fFactorFE < 1.0 ) || ( encodingPlan.fFactorPE < 1.0 )
    partialFourierText = sprintf(' (%dx%d)', ...
        encodingPlan.numFE, encodingPlan.numPE );
else
    partialFourierText = [];
end

% add text
textForImage.LBot{1} = ['FOV: ', ...
    num2str(1e3*imagePlane.fovX),'x',num2str(1e3*imagePlane.fovY),'mm'];
textForImage.LBot{2} = ['MATRIX: ', ...
    num2str(imagePlane.sizeX),'x',num2str(imagePlane.sizeY), ...
    partialFourierText];
textForImage.LBot{3} = sequenceText;
if ~isempty(irSenseText)
    textForImage.LBot{4} = irSenseText;
end
