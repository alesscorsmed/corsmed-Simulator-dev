function [sliceModel] = generateSlice( ...
    slicePlane, anatomicalModel, coilSystem, modelControl, dbgControl )
%
% SIMMODEL.HANDLING.GENERATESLICE
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
functionName = 'sliceModel:handling:generateSlice';
% check args
if (nargin < 4) || isempty(slicePlane) || isempty(anatomicalModel) ...
        || isempty(coilSystem) || isempty(modelControl)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 5) || isempty(dbgControl)
    dbgControl.mode = 0;
    dbgControl.file = [];
end
% info for debugging
if dbgControl.mode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(dbgControl.file,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end


%% decide which method is going to be used based on the requirements
if modelControl.isotropicGrid
    %% generation of an isotropic grid and interpolation into it
    [sliceModel] = simModel.handling.isotropicModelGeneration(...
        slicePlane, anatomicalModel, coilSystem, ...
        modelControl, dbgControl);
    
else
    %% non-isotropic: use of the original anatomical model grid
    [sliceModel] = simModel.handling.anisotropicModelGeneration(...
        slicePlane, anatomicalModel, coilSystem, ...
        modelControl, dbgControl );
end

%% TODO : clean / organize this
varyT1T2values = 0;
if varyT1T2values
    
    percT1variation = 0.1;
    percT2variation = 0.1;
    
    fprintf(1,'WARNING: Normal distribution of T1 and T2 values has been activated.\n')
    
    modelNew.tissueType         = 1:size(model.tissueType);
    modelNew.tissueType         = [1:size(model.tissueType)]';
    modelNew.tissueValues       = [zeros(size(model.tissueType,1),6)];
    modelNew.tissueValues(:,1)  = normrnd(model.tissueValues(model.tissueType,1),...
        percT1variation*model.tissueValues(model.tissueType,1));
    modelNew.tissueValues(:,2)  = normrnd(model.tissueValues(model.tissueType,2),...
        percT2variation*model.tissueValues(model.tissueType,2));
    modelNew.tissueValues(:,3)  = model.tissueValues(model.tissueType,3);
    modelNew.tissueValues(:,4)  = model.tissueValues(model.tissueType,4);
    modelNew.tissueValues(:,5)  = model.tissueValues(model.tissueType,5);
    modelNew.tissueValues(:,6)  = model.tissueValues(model.tissueType,6);
    
    model.tissueValues  = modelNew.tissueValues;
    model.tissueType    = modelNew.tissueType;
    
end

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

