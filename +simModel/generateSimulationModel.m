function [fovModel] = generateSimulationModel( ...
    fovDomain, anatomicalModel, coilSystem, modelControl, dbgControl )
%
% SIMMODEL.GENERATESIMULATIONMODEL
%
%     Function that generates the simulation model, by interpolating
%     data from the anatomical model and the coils into the disctretized
%     plane of the current slice
%
% INPUT
%
%
%
% OUTPUT
%   spinModel          updated spinModel struct with model
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'simModel:generateSimulationModel';
if (nargin < 3)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
if (nargin < 4) || isempty(modelControl)
    % initialize the simulation control w/ defaults
    [modelControl] = simModel.setup.initializeModelControl();
end
if (nargin < 5) || isempty(dbgControl)
    dbgControl.mode = 0;
    dbgControl.file = [];
else
    try
        %% TODO : fix expControl data
        dbgControl.mode = dbgControl.debugMode;
        dbgControl.file = dbgControl.debugFile;
    catch
        dbgControl.mode = 1;
        dbgControl.file = [];
    end
end

%% open debugging
if dbgControl.mode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(dbgControl.file,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n\n%s : start', functionName);
end

%% generate fovModel
fovModel.name       = fovDomain.name;

%% depending on 2D or 3D, generate the 'slice' Model
if fovDomain.is3D
    %% interpolate the 3D slab
    fovModel.is3D       = 1;
    fovModel.numSlices  = 1;
    fovModel.slice{1}.model = simModel.handling.generateSlice( ...
        fovDomain.slab.plane, anatomicalModel, coilSystem, ...
        modelControl, dbgControl);
else
    %% interpolate and generate simulation model for each slice
    fovModel.name       = fovDomain.name;
    fovModel.is3D       = 0;
    fovModel.numSlices  = fovDomain.numSlices;
    for ss = 1:fovModel.numSlices
        fovModel.slice{ss}.model = simModel.handling.generateSlice( ...
            fovDomain.slice{ss}.plane, anatomicalModel, coilSystem, ...
            modelControl, dbgControl);
    end
end
% frames and contrasts
if ~isfield(anatomicalModel,'numFrames')
    fovModel.numPhases  = 1;
else
    fovModel.numPhases	= anatomicalModel.numFrames;
end
fovModel.numContrasts   = 1;


%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : simulation model generated', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

