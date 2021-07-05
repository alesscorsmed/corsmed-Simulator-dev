function [fovDomain, fovModel] = generateSlices( ...
    acquisition, anatomicalModel, coilSystem, expControl)
%
% DOMAIN.GENERATESLICES
%
%     Function that handles the positions from Front End
%     and generates a structure with the simulation domain
%     including limit points, translations, rotations, ect...
%     for each slice and for the 3D domain
%     It also interpolates the model into the slices or slab of
%     the domain.
%
% INPUT
%   expControl         experiment control data, with commLocalDB
%   acquisition        structure with acquisition data
%
% OUTPUT
%   fovDomain          initialized struct with planes
%   fovModel           initialized struct with models
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'domain.generateSlices';
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

%% slice the domain to generate the planes
[fovDomain] = domain.generateSimulationDomain( acquisition, expControl );

%% depending on 2D or 3D, generate the Model
if fovDomain.is3D
    %% interpolate the 3D slab
    fovModel.name       = fovDomain.name;
    fovModel.is3D       = 1;
    fovModel.numSlices  = 1;
    fovModel.slice{1}.model = domain.generateSimulationModel( ...
        fovDomain.slab.plane, anatomicalModel, coilSystem, expControl);
else
    %% interpolate and generate simulation model for each slice
    fovModel.name       = fovDomain.name;
    fovModel.is3D       = 0;
    fovModel.numSlices  = fovDomain.numSlices;
    for ss = 1:fovModel.numSlices
        fovModel.slice{ss}.model = domain.generateSimulationModel( ...
            fovDomain.slice{ss}.plane, anatomicalModel, coilSystem, expControl);
    end
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for spinModel %s',...
        functionName, fovModel.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
        
