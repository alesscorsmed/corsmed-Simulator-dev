function [BOrient,TOrient,LOrient,ROrient, planeType, imageType] = ...
    planeOrientation( pLBot, pRBot, pLTop, pRTop)
%
% DOMAIN.PLANEHANDLING.PLANEORIENTATION
%
%     Gets the coordinates of the plane and determines which is the 
%     actual orientation of the plane (position of the points)
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

pLBot = reshape(pLBot,3,1);
pRBot = reshape(pRBot,3,1);
pLTop = reshape(pLTop,3,1);
pRTop = reshape(pRTop,3,1);

%% unitatary:
I = eye(3);

%% get the normal of the plane
% use Left Bottom as Reference
vecBottom  = pRBot - pLBot;
vecLeft    = pLTop - pLBot;
normalVec  = cross(vecBottom, vecLeft);

%% get the projections on each unitary vector:
normalProj = I*normalVec;
% find closest axis
[~, normalIdx] = max(abs(normalProj));

%% define the type of plane and image
% depending on the dominant direction, the type of image is
switch normalIdx
    case 1
        imageType = 'Sagittal';
    case 2
        imageType = 'Coronal';
    case 3
        imageType = 'Transversal';
end
% depending on the direction of the normal, type of plane
switch nnz(normalProj)
    case 3
        planeType = 'double-oblique';
    case 2
        planeType = 'oblique';
    case 1
        planeType = 'in-plane';
end

%% apply rotation to take to the closest plane
% unitary vector axis
normalAxis = sign(normalProj(normalIdx))*I(:,normalIdx);
% find rotation to move the normal to that axis
rot1Vec = vrrotvec(normalProj,normalAxis);
rot1Mat = vrrotvec2mat(rot1Vec);
% rotate the points of the plane
rot1LBot = rot1Mat*pLBot;
rot1RBot = rot1Mat*pRBot;
rot1LTop = rot1Mat*pLTop;
rot1RTop = rot1Mat*pRTop;

%% find second rotation to bring the (new) bottom vector to the closest axis
planeVec  = rot1RBot - rot1LBot;
planeProj = I*planeVec;
% find closest axis
[~, planeIdx] = max(abs(planeProj));
% depending on the dominant direction, bottom axis is
switch planeIdx
    case 1
        bottomDir = 'LR';
    case 2
        bottomDir = 'AP';
    case 3
        bottomDir = 'HF';
end
%% apply rotation to take to the closest axis
% unitary vector axis
planeAxis = sign(planeProj(planeIdx))*I(:,planeIdx);
% find rotation to move the normal to that axis
rot2Vec = vrrotvec(planeProj,planeAxis);
rot2Mat = vrrotvec2mat(rot2Vec);
% rotate the points of the plane
rot2LBot = rot2Mat*rot1LBot;
rot2RBot = rot2Mat*rot1RBot;
rot2LTop = rot2Mat*rot1LTop;
rot2RTop = rot2Mat*rot1RTop;

% ORICOORD    = [pLBot, pRBot, pRTop, pLTop, pLBot].';
% ORIVEC      = [pLBot, pRBot].';
% ROT1COORD   = [rot1LBot, rot1RBot, rot1RTop, rot1LTop, rot1LBot].';
% ROT1VEC     = [rot1LBot, rot1RBot].';
% ROT2COORD   = [rot2LBot, rot2RBot, rot2RTop, rot2LTop, rot2LBot].';
% ROT2VEC     = [rot2LBot, rot2RBot].';
%     
% plot3(ORICOORD(:,1), ORICOORD(:,2), ORICOORD(:,3), 'r--');
% hold on;
% plot3(ORIVEC(:,1), ORIVEC(:,2), ORIVEC(:,3), 'r-');
% plot3(pLBot(1), pLBot(2), pLBot(3), 'ro');
% hold on;
% plot3(ROT1COORD(:,1), ROT1COORD(:,2), ROT1COORD(:,3), 'g--');
% plot3(ROT1VEC(:,1), ROT1VEC(:,2), ROT1VEC(:,3), 'g-');
% plot3(rot1LBot(1), rot1LBot(2), rot1LBot(3), 'go');
% hold on;
% plot3(ROT2COORD(:,1), ROT2COORD(:,2), ROT2COORD(:,3), 'b--');
% plot3(ROT2VEC(:,1), ROT2VEC(:,2), ROT2VEC(:,3), 'b-');
% plot3(rot2LBot(1), rot2LBot(2), rot2LBot(3), 'bo');
% grid on;

%% define now the position of each point
% see what is the directtion of the bottom line
vecB = rot2RBot - rot2LBot;
[~, idxB] = max(abs(vecB));
switch idxB   
    case 1 % bottom line runs along X (LR)
        if vecB(idxB) < 0
            LOrient = 'L';
            ROrient = 'R';
        else
            LOrient = 'R';
            ROrient = 'L';
        end
    case 2 % bottom line runs along Y (PA)
        if vecB(idxB) < 0
            LOrient = 'P';
            ROrient = 'A';
        else
            LOrient = 'A';
            ROrient = 'P';
        end
    case 3 % bottom line runs along Z (HF)
        if vecB(idxB) < 0
            LOrient = 'H';
            ROrient = 'F';
        else
            LOrient = 'F';
            ROrient = 'H';
        end
end
% see what is the direction of the left line
vecL = rot2LTop - rot2LBot;
[~, idxL] = max(abs(vecL));
switch idxL   
    case 1 % line runs along X (LR)
        if vecL(idxL) < 0
            BOrient = 'L';
            TOrient = 'R';
        else
            BOrient = 'R';
            TOrient = 'L';
        end
    case 2 % line runs along Y (PA)
        if vecL(idxL) < 0
            BOrient = 'P';
            TOrient = 'A';
        else
            BOrient = 'A';
            TOrient = 'P';
        end
    case 3 % line runs along Z (FH)
        if vecL(idxL) < 0
            BOrient = 'H';
            TOrient = 'F';
        else
            BOrient = 'F';
            TOrient = 'H';
        end
end
% 
% fid = 1;
% fprintf(fid, '\n NEW');
% fprintf(fid, '\n  Top    Left       [%.3f, %.3f, %.3f ]', pLTop);
% fprintf(fid, '\n  Top    Right      [%.3f, %.3f, %.3f ]', pRTop);
% fprintf(fid, '\n  Bottom Right      [%.3f, %.3f, %.3f ]', pRBot);
% fprintf(fid, '\n  Bottom Left       [%.3f, %.3f, %.3f ]', pLBot);
% fprintf(fid, '\n  Plane Type        %s', planeType);
% fprintf(fid, '\n  Image Type        %s', imageType);
% fprintf(fid, '\n                        %s', TOrient);
% fprintf(fid, '\n  Orientation        %s     %s', LOrient, ROrient);
% fprintf(fid, '\n                        %s', BOrient);
% fprintf(fid, '\n ');

