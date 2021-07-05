function [plane] = generatePlane(plane)
%
% generates a centered plane from slice to simulate
%

% handle the plane to get the positions
[p1OO, p2FE, p3PE, LBotNew, RBotNew, LTopNew, RTopNew] = ...
    domain.planeHandling.planePosition( ...
    plane.LTopFrontEnd, plane.RTopFrontEnd, ...
    plane.LBotFrontEnd, plane.RBotFrontEnd, ...
    plane.LTop, plane.RTop, ...
    plane.LBot, plane.RBot );

% Calculate the rotation matrices based on the 3 points
[rotMatX, rotMatY, rotMatZ, p2Final, p3Final] = ...
    domain.planeHandling.planeRotation(p1OO, p2FE, p3PE);

% The slice is also translated now (the ref point for the rotations is p1!)
transX = p2Final(1)/2;
transY = p3Final(2)/2;
% generate the plane, given by 3 coorners
planeCoords = [ zeros(3,1), p2Final, p3Final].';
% move the center of the plane to (0,0,0)
planeCoords(:,1) = planeCoords(:,1) - transX;
planeCoords(:,2) = planeCoords(:,2) - transY;
planeCoords(:,3) = 0;

% orientations
try
    [BOrient, TOrient, LOrient, ROrient, planeType, imageType] = ...
        domain.planeHandling.planeOrientation(LBotNew, RBotNew, LTopNew, RTopNew);
catch
    BOrient = [];
    TOrient = [];
    LOrient = [];
    ROrient = [];
    planeType = [];
    imageType = [];
end

% assign positions
plane.p1OO = p1OO;
plane.p2FE = p2FE;
plane.p3PE = p3PE;
%
plane.p2Final = p2Final;
plane.p3Final = p3Final;
% rots and trans
plane.rotMatX = rotMatX;
plane.rotMatY = rotMatY;
plane.rotMatZ = rotMatZ;
%
plane.transX = transX;
plane.transY = transY;
% coords and orientations
plane.planeCoords = planeCoords;
%
plane.BOrient = BOrient;
plane.TOrient = TOrient;
plane.LOrient = LOrient;
plane.ROrient = ROrient;
plane.Type    = planeType;
plane.imType  = imageType;
%
plane.LBotNew = LBotNew;
plane.RBotNew = RBotNew;
plane.LTopNew = LTopNew;
plane.RTopNew = RTopNew;