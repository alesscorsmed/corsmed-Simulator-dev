function [spinModel] = generateSimulationDomain2D(spinModel,acqData,expControl)
%
% DOMAIN.GENERATESIMULATIONDOMAIN2D
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
%   spinModel          initialized spinModel struct
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'domain.generateSimulationDomain2D';
if (nargin < 2)
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

%% get all the plane points
pointsAll           = acqData.pointsAll;
pointsAllFrontend   = acqData.pointsAllFrontend;
assert(spinModel.numSlices==length(pointsAll),...
    sprintf('\nERROR: %s : number of slices do not match plane points \n',functionName));

%% get 3D thickness
spinModel.is3D = 0;
numSlices      = acqData.numSlices;
sliceThickness = acqData.sliceThickness;
sliceGap       = acqData.sliceGap;
slabThickness  = numSlices*sliceThickness + (numSlices-1)*sliceGap;

%% populate the spinModel plane with slice data
for ss = 1:numSlices
    % extract point data
    points = pointsAll{ss};
    points = strsplit(points,';');
    % front end points
    pointsFrontEnd  = pointsAllFrontend{ss};
    pointsFrontEnd  = strsplit(pointsFrontEnd,';');
    % assig points
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
    % assign the center position and slice thickness
    plane.thickness = sliceThickness;
    plane.zPosition = (sliceThickness-slabThickness)/2 ...
        + (ss-1)*(sliceThickness + sliceGap); % start from bottom
    % assign to slice and spinModel
    spinModel.slice{ss}.plane = plane;
end

% Generate the center plane
if mod(spinModel.numSlices,2)
    % If odd number of slices, take the points of the middle slice
    ss = round(spinModel.numSlices/2);
    plane = spinModel.slice{ss}.plane;
    % assign the slice thickness
    plane.zPosition = 0.0;
    plane.thickness = slabThickness;
    % assign to slice and spinModel
    spinModel.center.plane = plane;
else
    % average center slices
    s1 = floor(spinModel.numSlices/2);
    s2 = s1 + 1;
    % assig points
    plane.LTop = (spinModel.slice{s1}.plane.LTop + ...
        spinModel.slice{s2}.plane.LTop)/2;
    plane.RTop = (spinModel.slice{s1}.plane.RTop + ...
        spinModel.slice{s2}.plane.RTop)/2;
    plane.LBot = (spinModel.slice{s1}.plane.LBot + ...
        spinModel.slice{s2}.plane.LBot)/2;
    plane.RBot = (spinModel.slice{s1}.plane.RBot + ...
        spinModel.slice{s2}.plane.RBot)/2;   
    % get the points from the Front End
    plane.LTopFrontEnd = (spinModel.slice{s1}.plane.LTopFrontEnd + ...
        spinModel.slice{s2}.plane.LTopFrontEnd)/2;
    plane.RTopFrontEnd = (spinModel.slice{s1}.plane.RTopFrontEnd + ...
        spinModel.slice{s2}.plane.RTopFrontEnd)/2;
    plane.LBotFrontEnd = (spinModel.slice{s1}.plane.LBotFrontEnd + ...
        spinModel.slice{s2}.plane.LBotFrontEnd)/2;
    plane.RBotFrontEnd = (spinModel.slice{s1}.plane.RBotFrontEnd + ...
        spinModel.slice{s2}.plane.RBotFrontEnd)/2;  
    % compute the transformations and final plane coordinates
    [plane] = domain.planeHandling.generatePlane(plane);
    % assign the slice thickness
    plane.zPosition = 0.0;
    plane.thickness = slabThickness;
    % assign to slice and spinModel
    spinModel.center.plane = plane;
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done ', functionName );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
        
