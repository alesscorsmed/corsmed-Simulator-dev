function [spinModel] = generateSimulationDomain3D(spinModel,acqData,expControl)
%
% DOMAIN.GENERATESIMULATIONDOMAIN3D
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
functionName = 'domain.generateSimulationDomain3D';
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
assert(acqData.numSlices==size(pointsAll,2),...
    sprintf('\nERROR: %s : number of slices do not match plane points \n',functionName));

%% get 3D thickness
spinModel.is3D = 1;
numSlices      = acqData.numSlices;
sliceThickness = acqData.sliceThickness;
sliceGap       = acqData.sliceGap;
slabThickness  = numSlices*sliceThickness + (numSlices-1)*sliceGap;

%% populate the spinModel plane with slice data
% Generate the center plane
if mod(numSlices,2)
    % If odd number of slices, take the points of the middle slice
    ss = round(numSlices/2);
    % extract point data
    points = pointsAll{1,ss};
    points = strsplit(points,';');
    % front end points
    pointsFrontEnd  = pointsAllFrontend{1,ss};
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
    % assign the slice thickness
    plane.zPosition = 0.0;
    plane.thickness = slabThickness;
    % assign to slice and spinModel
    spinModel.slice{1}.plane = plane;
    spinModel.center.plane = plane;
else
    % average center slices
    s1 = floor(numSlices/2);
    s2 = s1 + 1;
    % extract points
    pointsS1 = pointsAll{1,s1};
    pointsS1 = strsplit(pointsS1,';');
    pointsS2 = pointsAll{1,s2};
    pointsS2 = strsplit(pointsS2,';');
    % front end points
    pointsFrontEndS1  = pointsAllFrontend{1,s1};
    pointsFrontEndS1  = strsplit(pointsFrontEndS1,';');
    pointsFrontEndS2  = pointsAllFrontend{1,s2};
    pointsFrontEndS2  = strsplit(pointsFrontEndS2,';');
    % assig averga points
    plane.LTop = ( str2num(pointsS1{1,1}) + str2num(pointsS2{1,1}) )/2;
    plane.RTop = ( str2num(pointsS1{1,2}) + str2num(pointsS2{1,2}) )/2;
    plane.LBot = ( str2num(pointsS1{1,3}) + str2num(pointsS2{1,3}) )/2;
    plane.RBot = ( str2num(pointsS1{1,4}) + str2num(pointsS2{1,4}) )/2;   
    % get the points from the Front End
    plane.LTopFrontEnd = ( str2num(pointsFrontEndS1{1,1}) ...
        + str2num(pointsFrontEndS2{1,1}) )/2;
    plane.RTopFrontEnd = ( str2num(pointsFrontEndS1{1,2}) ...
        + str2num(pointsFrontEndS2{1,2}) )/2;
    plane.LBotFrontEnd = ( str2num(pointsFrontEndS1{1,3}) ...
        + str2num(pointsFrontEndS2{1,3}) )/2;
    plane.RBotFrontEnd = ( str2num(pointsFrontEndS1{1,4}) ...
        + str2num(pointsFrontEndS2{1,4}) )/2; 
    % compute the transformations and final plane coordinates
    [plane] = domain.planeHandling.generatePlane(plane);
    % assign the slice thickness
    plane.zPosition = 0.0;
    plane.thickness = slabThickness;
    % assign to slice and spinModel
    spinModel.slice{1}.plane = plane;
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
        
