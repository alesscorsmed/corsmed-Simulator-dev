function [reconData] = generateSlices( ...
    fovDomain, anatomicalModel, coilSystem, ...
    encodingData, noiseData, expControl )
%
% RECONSTRUCTOR.GENERATESLICES
%
%     Initializes the reconData and populates the 
%     slices with data
%
% INPUT
%   fovDomain          struct with domain planes
%   anatomicalModel    anatomical model
%   expControl         experiment control data
%
% OUTPUT
%   reconData         initialized struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructor.generateSlices';
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

%% assign data to reconData
reconData.name	= fovDomain.name;
reconData.is3D  = fovDomain.is3D;
% % image sizes
% reconData.imSizeX = encodingData.plan.imSizeX;
% reconData.imSizeY = encodingData.plan.imSizeY;
% reconData.imSizeZ = encodingData.plan.imSizeZ;
% % FOV of image
% reconData.imFovX = encodingData.plan.imFovX;
% reconData.imFovY = encodingData.plan.imFovY;
% reconData.imFovZ = encodingData.plan.imFovZ;
% k-space sizes
reconData.numX = 0;
reconData.numY = 0;
reconData.numZ = 0;
reconData.numC = 0;

%% add encoding data
reconData.encoding      = encodingData;

%% and noise data
noiseData.voxelVolume           = prod(expControl.model.gridStep);
noiseData.defaultVoxelVolume    = prod([0.001,0.001,0.001]);
% Calculate noiseAmp so as give similar noise level regardless the 
% selected grid size
noiseData.noiseAmp              = (noiseData.noiseSD*noiseData.defaultVoxelVolume);
reconData.noise                 = noiseData;

%% prepare frames
reconData.numFrames     = encodingData.plan.numFrames;
reconData.numCtFr       = encodingData.plan.numContrasts; % contrast frames
reconData.numPhFr       = encodingData.plan.numPhases; % phase frames
ff = 0;
for cf = 1:reconData.numCtFr
    for pf = 1:reconData.numPhFr
        ff = ff+1; % increase frame counter
        %% depending on 2D or 3D, generate the Model
        if fovDomain.is3D
            %% interpolate the 3D slab
            reconData.numSlices  = 1;
            % interpolate tissue mask and Coil sens in the image coordinates
            [tissueMask, Cx, Cy] = reconstructor.tools.generateMaps( ...
                fovDomain.slab.plane, encodingData.plan, ...
                anatomicalModel, coilSystem, expControl );
            % assign
            reconData.slice{1}.sens               = Cx - 1j*Cy; % sensitivities
            reconData.slice{1}.mask               = tissueMask; % tissue mask
            reconData.slice{1}.frame{ff}.raw      = []; % struct with time signal data, assembled from multiple simSignal
            reconData.slice{1}.frame{ff}.kSpace   = []; % FE x PE x SE x COILS x CONTRASTS
            reconData.slice{1}.frame{ff}.iSpace   = []; % FE x PE x SE x   1   x CONTRASTS
            reconData.slice{1}.frame{ff}.phIdx    = pf;
            reconData.slice{1}.frame{ff}.ctIdx    = cf;
        else
            %% interpolate and generate simulation model for each slice
            reconData.numSlices = fovDomain.numSlices;
            for ss = 1:reconData.numSlices
                % interpolate tissue mask and Coil sens in the image coordinates
                [tissueMask, Cx, Cy] = reconstructor.tools.generateMaps( ...
                    fovDomain.slice{ss}.plane, encodingData.plan, ...
                    anatomicalModel, coilSystem, expControl );
                % assign
                reconData.slice{ss}.sens               = Cx - 1j*Cy; % sensitivities
                reconData.slice{ss}.mask               = tissueMask; % tissue mask
                reconData.slice{ss}.frame{ff}.raw      = []; % struct with time signal data, assembled from multiple simSignal
                reconData.slice{ss}.frame{ff}.kSpace   = [];% FE x PE x SE x COILS x CONTRASTS
                reconData.slice{ss}.frame{ff}.iSpace   = [];% FE x PE x SE x   1   x CONTRASTS
                reconData.slice{ss}.frame{ff}.phIdx    = pf;
                reconData.slice{ss}.frame{ff}.ctIdx    = cf;
            end
        end

    end
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for model %s',...
        functionName, reconData.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end