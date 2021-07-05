function [imageData] = generateSlices( ...
    fovDomain, pulseSequence, encodingPlan, expControl)
%
% IMAGE.GENERATESLICES
%
%     Initializes the imageData and populates the 
%     slices with coordinates and text
%
% INPUT
%   fovDomain          struct with domain planes
%   anatomicalModel    anatomical model
%   expControl         experiment control data
%
% OUTPUT
%   imageData           initialized struct with models
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'image.generateSlices';
if (nargin < 3)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n\n%s : start', functionName);
end

%% initialize imageData
imageData.name          = fovDomain.name;
imageData.bodyPartName  = fovDomain.bodyPartName;
imageData.is3D          = fovDomain.is3D;
imageData.numSlices     = fovDomain.numSlices;
imageData.numContrasts  = encodingPlan.numCE;
imageData.b0map         = expControl.b0map;
imageData.numFrames     = encodingPlan.numFrames;
imageData.numCtFr       = encodingPlan.numContrasts; % contrast frames
imageData.numPhFr       = encodingPlan.numPhases; % phase frames

%% define number of contrast and post-processing method based on sequence
[imageData] = image.dicom.defineContrastProcess( ...
    imageData, pulseSequence.familyName );

%% prepare the frames
ff = 0;
for cf = 1:imageData.numCtFr
    for pf = 1:imageData.numPhFr
        ff = ff+1; % increase frame counter
        %% prepare the slices
        for ss = 1:imageData.numSlices
            
            %% generate the image planes
            imagePlane = image.plane.generateImagePlane( ...
                fovDomain.slice{ss}.plane, expControl.outerFOVratio);
            % add image size in pixels
            imagePlane.sizeX = encodingPlan.imSizeX;
            imagePlane.sizeY = encodingPlan.imSizeY;
            % assign to image slice
            imageSlice.plane = imagePlane;
            
            %% for each slice, allocate contrasts
            for cc = 1:imageData.numContrasts
                imageSlice.frame{ff}.phIdx = pf;
                imageSlice.frame{ff}.ctIdx = cf;
                imageSlice.frame{ff}.contrast{cc}.uniqueID    = [];
                imageSlice.frame{ff}.contrast{cc}.image       = [];
                imageSlice.frame{ff}.contrast{cc}.mask        = [];
                % generate the text
                imageSlice.frame{ff}.contrast{cc}.imageText = image.dicom.prepareViewerText( ...
                    imagePlane, imageData, pulseSequence, ...
                    encodingPlan, expControl, ss, cc );
                
            end
            %% assign
            imageData.slice{ss} = imageSlice;
        end
    end
end

%% Transfer in imageData pulse-sequence-specific data
if strcmpi(pulseSequence.familyName,'molli')
    imageData.sequence.MOLLI = pulseSequence.MOLLI;
end


%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for model %s',...
        functionName, imageData.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
