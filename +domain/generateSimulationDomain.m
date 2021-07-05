function [fovDomain] = generateSimulationDomain(acquisition,expControl)
%
% DOMAIN.GENERATESIMULATIONDOMAIN
%
%     Function that handles the positions from Front End
%     and generates a structure with the simulation domain
%     including limit points, translations, rotations, ect...
%     for each slice and for the 3D domain
%
% INPUT
%   acquisition        structure with acquisition data
%   expControl         experiment control data, with commLocalDB
%
% OUTPUT
%   fovDomain          initialized fovDomain struct
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'domain.generateSimulationDomain';
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
    fprintf(fid, '\n%s : start', functionName);
end

%% initialize the domain with the number of slices and number of phases
acqData     = acquisition.data;
numSlices   = acqData.numSlices;

%% get all the plane points
pointsAll           = acqData.pointsAll;
pointsAllFrontend   = acqData.pointsAllFrontend;
assert(numSlices==length(pointsAll),...
    sprintf('\nERROR: %s : number of slices do not match plane points \n',functionName));

%% get 3D thickness
sliceThickness   	= acqData.sliceThickness;
sliceGap          	= acqData.sliceGap;
slabThickness     	= numSlices*sliceThickness + (numSlices-1)*sliceGap;
fovDomain.is3D      = acqData.is3D;
fovDomain.numSlices = numSlices;

%% populate the fovDomain plane with slice data
for ss = 1:numSlices
    
    %% extract coordinates data
    points          = pointsAll{ss};
    points          = strsplit(points,';');
    pointsFrontEnd  = pointsAllFrontend{ss};
    pointsFrontEnd  = strsplit(pointsFrontEnd,';');
    
    %% generate plane
    % plane FOV
    plane.fovX = acqData.fovFE;
    plane.fovY = acqData.fovPE;
    plane.foldoverDir = acqData.foldoverDir;
    % assign points
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
    % assign the (rotated) reference position and slice thickness
    plane.refRot = plane.p1OO*(plane.rotMatX*plane.rotMatY*plane.rotMatZ).';
    plane.thickness = sliceThickness;
    % assign to slice in fovDomain
    fovDomain.slice{ss}.plane = plane;
    
end

%% Generate the center plane
if mod(fovDomain.numSlices,2)
    % If odd number of slices, take the points of the middle slice
    ss = round(fovDomain.numSlices/2);
    plane = fovDomain.slice{ss}.plane;
    % plane FOV
    plane.fovX = acqData.fovFE;
    plane.fovY = acqData.fovPE;
    % assign the slice thickness
    plane.zPosition = 0.0;
    plane.thickness = slabThickness;
    % assign to slice and fovDomain
    fovDomain.slab.plane = plane;
else
    % average center slices
    s1 = floor(fovDomain.numSlices/2);
    s2 = s1 + 1;
    % assig points
    plane.LTop = (fovDomain.slice{s1}.plane.LTop + ...
        fovDomain.slice{s2}.plane.LTop)/2;
    plane.RTop = (fovDomain.slice{s1}.plane.RTop + ...
        fovDomain.slice{s2}.plane.RTop)/2;
    plane.LBot = (fovDomain.slice{s1}.plane.LBot + ...
        fovDomain.slice{s2}.plane.LBot)/2;
    plane.RBot = (fovDomain.slice{s1}.plane.RBot + ...
        fovDomain.slice{s2}.plane.RBot)/2;   
    % get the points from the Front End
    plane.LTopFrontEnd = (fovDomain.slice{s1}.plane.LTopFrontEnd + ...
        fovDomain.slice{s2}.plane.LTopFrontEnd)/2;
    plane.RTopFrontEnd = (fovDomain.slice{s1}.plane.RTopFrontEnd + ...
        fovDomain.slice{s2}.plane.RTopFrontEnd)/2;
    plane.LBotFrontEnd = (fovDomain.slice{s1}.plane.LBotFrontEnd + ...
        fovDomain.slice{s2}.plane.LBotFrontEnd)/2;
    plane.RBotFrontEnd = (fovDomain.slice{s1}.plane.RBotFrontEnd + ...
        fovDomain.slice{s2}.plane.RBotFrontEnd)/2;  
    % compute the transformations and final plane coordinates
    [plane] = domain.planeHandling.generatePlane(plane);
    % plane FOV
    plane.fovX = acqData.fovFE;
    plane.fovY = acqData.fovPE;
    % assign the (rotated) reference position andslice thickness 
    plane.refRot = plane.p1OO*(plane.rotMatX*plane.rotMatY*plane.rotMatZ).';
    plane.zPosition = 0.0;
    plane.thickness = slabThickness;
    % assign to slice and fovDomain
    fovDomain.slab.plane = plane;
end

%% find the actual relative z position of each slice w.r.t. center
fovDomain.slab.zSliceCoord = zeros(numSlices,1);
for ss = 1:numSlices
    fovDomain.slice{ss}.plane.zPosition = ...
        fovDomain.slice{ss}.plane.refRot(3) - fovDomain.slab.plane.refRot(3);
    %fovDomain.slice{ss}.plane.zPosition2 = (sliceThickness-slabThickness)/2 ...
    %    + (ss-1)*(sliceThickness + sliceGap); % start from bottom
    fovDomain.slab.zSliceCoord(ss) = fovDomain.slice{ss}.plane.zPosition;
end

%% name
switch expControl.courseID
    case 4
        fovDomain.name          = sprintf('XCAT-ExpID%d',expControl.experimentID);
        fovDomain.bodyPartName  = 'Body';
    case 6
        fovDomain.name = sprintf('MIDA-ExpID%d',expControl.experimentID);
        fovDomain.bodyPartName  = 'Head';
    case 8
        fovDomain.name = sprintf('VOXELMAN-ExpID%d',expControl.experimentID);
        fovDomain.bodyPartName  = 'Body';
    otherwise
        fovDomain.name = sprintf('UNKNOWN-ExpID%d',expControl.experimentID);
        fovDomain.bodyPartName  = 'Body';
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for model %s',...
        functionName, fovDomain.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    plane = fovDomain.slab.plane;
    fprintf(fid, '\n  Number of Slices  %d', fovDomain.numSlices);
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
    if ~isempty(plane.imType)
        fprintf(fid, '\n  Image Type        %s', plane.imType);
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
        
