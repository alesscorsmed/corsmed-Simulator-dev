function [model] = generateSimulationModel( ...
    plane, anatomicalModel, coilSystem, expControl )
%
% DOMAIN.GENERATESIMULATIONMODEL
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
functionName = 'domain.generateSimulationModel';
if (nargin < 4)
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

%% decide which method is going to be used based on the requirements
if expControl.model.isotropicGrid
    %% generation of an isotropic grid and interpolation into it
    [model] = domain.modelHandling.isotropicModelGeneration(...
        plane, anatomicalModel, coilSystem, expControl);
    
else
    %% non-isotropic: use of the original anatomical model grid
    [model] = domain.modelHandling.anisotropicModelGeneration(...
        plane, anatomicalModel, coilSystem, expControl );
end

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
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : simulation model generated', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

