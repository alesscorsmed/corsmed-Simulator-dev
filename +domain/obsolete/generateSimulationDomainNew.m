function [fovDomain] = generateSimulationDomain( ...
    acquisition, expControl)
%
% DOMAIN.GENERATESIMULATIONDOMAIN
%
%     Function that handles the positions from Front End
%     and generates a structure with the simulation domain
%     including limit points, translations, rotations, ect...
%     for each slice and for the 3D domain
%
% INPUT
%   expControl         experiment control data, with commLocalDB
%   acquisition        structure with acquisition data
%
% OUTPUT
%   simulationDomain   initialized struct with slices
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'domain.generateSimu';
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

%% initialize the spin model with the number of slices
if acquisition.data.is3D
    [spinModel] = data.spinModel.initialize(1);
    [spinModel] = domain.generateSimulationDomain3D( spinModel, ...
        acquisition.data, expControl);
else
    [spinModel] = data.spinModel.initialize(acquisition.data.numSlices);
    [spinModel] = domain.generateSimulationDomain2D( spinModel, ...
        acquisition.data, expControl);
end

%% interpolate and generate simulation model
for ss = 1:spinModel.numSlices
    spinModel.slice{ss}.model = domain.generateSimulationModel( ...
        spinModel.slice{ss}.plane, anatomicalModel, coilSystem, expControl);
end

%% name
switch expControl.courseID
    case 4
        spinModel.name          = sprintf('XCAT-ExpID%d',expControl.experimentID);
        spinModel.bodyPartName  = 'Body';
    case 6
        spinModel.name = sprintf('MIDA-ExpID%d',expControl.experimentID);
        spinModel.bodyPartName  = 'Head';
    case 8
        spinModel.name = sprintf('VOXELMAN-ExpID%d',expControl.experimentID);
        spinModel.bodyPartName  = 'Body';
    otherwise
        spinModel.name = sprintf('UNKNOWN-ExpID%d',expControl.experimentID);
        spinModel.bodyPartName  = 'Body';
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for spinModel %s',...
        functionName, spinModel.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    plane = spinModel.center.plane;
    fprintf(fid, '\n  Number of Slices  %d', spinModel.numSlices);
    fprintf(fid, '\n  Volumne Thickness %.3fmm', plane.thickness*1e3);
    fprintf(fid, '\n  Top    Left       [%.3f, %.3f, %.3f ]', plane.LTop);
    fprintf(fid, '\n  Top    Right      [%.3f, %.3f, %.3f ]', plane.RTop);
    fprintf(fid, '\n  Bottom Right      [%.3f, %.3f, %.3f ]', plane.RBot);
    fprintf(fid, '\n  Bottom Left       [%.3f, %.3f, %.3f ]', plane.LBot);
    fprintf(fid, '\n  Top    Left  (FE) [%.3f, %.3f, %.3f ]', plane.LTopFrontEnd);
    fprintf(fid, '\n  Top    Right (FE) [%.3f, %.3f, %.3f ]', plane.RTopFrontEnd);
    fprintf(fid, '\n  Bottom Right (FE) [%.3f, %.3f, %.3f ]', plane.RBotFrontEnd);
    fprintf(fid, '\n  Bottom Left  (FE) [%.3f, %.3f, %.3f ]', plane.LBotFrontEnd);
    if ~isempty(plane.Type)
        fprintf(fid, '\n  Plane Type        %s', plane.Type);
    end
    if ~isempty(plane.BOrient) && ~isempty(plane.TOrient) ...
            && ~isempty(plane.LOrient) && ~isempty(plane.ROrient)
        fprintf(fid, '\n                        %s', plane.TOrient);
        fprintf(fid, '\n  Orientation        %s     %s', plane.LOrient, plane.ROrient);
        fprintf(fid, '\n                        %s', plane.BOrient);
    end
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
        
