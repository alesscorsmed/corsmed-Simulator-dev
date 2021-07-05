function [imagePlane] = generateImagePlane( plane, outerFOVratio )
%
% IMAGE.PLANE.GENERATEIMAGEPLANE
%
%     Corrects the orientation of the image and generates coordinates
%
% INPUT
%   
%
% OUTPUT
%   
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'image.plane.generateImagePlane';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% extract FOV plane info
LTopFOV     = plane.LTopNew;
RTopFOV     = plane.RTopNew;
LBotFOV     = plane.LBotNew;
RBotFOV     = plane.RBotNew;
TOrientFOV 	= plane.TOrient;
BOrientFOV  = plane.BOrient;
ROrientFOV 	= plane.ROrient;
LOrientFOV 	= plane.LOrient;
fovX        = plane.fovX;
fovY        = plane.fovY;

%% get the correct orientation
dimensions = strcat(BOrientFOV,TOrientFOV,ROrientFOV,LOrientFOV);
dimensions = sort(dimensions);

switch dimensions
    case 'AFHP'
        TOrientImage = 'H';
        BOrientImage = 'F';
        LOrientImage = 'A';
        ROrientImage = 'P';
        
        if(TOrientFOV=='A')
            if(ROrientFOV=='H')
                LTopImage = RTopFOV;
                RTopImage = RBotFOV;
                LBotImage = LTopFOV;
                RBotImage = LBotFOV;
                imageFlip = 0;
                imageRot  = 90;
            else
                % mirror + rotate counterclock
                LTopImage = LTopFOV;
                RTopImage = LBotFOV;
                LBotImage = RTopFOV;
                RBotImage = RBotFOV;
                imageFlip = 2;
                imageRot  = 90;
            end
            lengthX = fovY;
            lengthY = fovX;
        elseif(TOrientFOV=='F')
            if(ROrientFOV=='A')
                % rotate counterclock 2
                LTopImage = RBotFOV;
                RTopImage = LBotFOV;
                LBotImage = RTopFOV;
                RBotImage = LTopFOV;
                imageFlip = 0;
                imageRot  = 180;
            else
                % mirror + rotate counterclock 2
                LTopImage = LBotFOV;
                RTopImage = RBotFOV;
                LBotImage = LTopFOV;
                RBotImage = RTopFOV;
                imageFlip = 2;
                imageRot  = 180;
            end
            lengthX = fovX;
            lengthY = fovY;
        elseif(TOrientFOV=='P')
            if(ROrientFOV=='F')
                % rotate clock wise 1
                LTopImage = LBotFOV;
                RTopImage = LTopFOV;
                LBotImage = RBotFOV;
                RBotImage = RTopFOV;
                imageFlip = 0;
                imageRot  = -90;
            else
                % mirror + clock wise 1
                LTopImage = RBotFOV;
                RTopImage = RTopFOV;
                LBotImage = LBotFOV;
                RBotImage = LTopFOV;
                imageFlip = 2;
                imageRot  = -90;
            end
            lengthX = fovY;
            lengthY = fovX;
        else
            if(ROrientFOV=='P')
                % Proper Orientation
                LTopImage = LTopFOV;
                RTopImage = RTopFOV;
                LBotImage = LBotFOV;
                RBotImage = RBotFOV;
                imageFlip = 0;
                imageRot  = 0;
            else
                % mirror
                LTopImage = RTopFOV;
                RTopImage = LTopFOV;
                LBotImage = RBotFOV;
                RBotImage = LBotFOV;
                imageFlip = 2;
                imageRot  = 0;
            end
            lengthX = fovX;
            lengthY = fovY;
        end
        
    case 'ALPR'
        TOrientImage = 'A';
        BOrientImage = 'P';
        LOrientImage = 'R';
        ROrientImage = 'L';
        
        if(TOrientFOV=='A')
            if(ROrientFOV=='L')
                % Proper Orientation
                LTopImage = LTopFOV;
                RTopImage = RTopFOV;
                LBotImage = LBotFOV;
                RBotImage = RBotFOV;
                imageFlip = 0;
                imageRot  = 0;
            else
                % mirror
                LTopImage = RTopFOV;
                RTopImage = LTopFOV;
                LBotImage = RBotFOV;
                RBotImage = LBotFOV;
                imageFlip = 2;
                imageRot  = 0;
            end
            lengthX = fovX;
            lengthY = fovY;
        elseif(TOrientFOV=='L')
            if(ROrientFOV=='P')
                % rotate clockwise 1
                LTopImage = LBotFOV;
                RTopImage = LTopFOV;
                LBotImage = RBotFOV;
                RBotImage = RTopFOV;
                imageFlip = 0;
                imageRot  = -90;
            else
                % mirror + clockwise 1
                LTopImage = RBotFOV;
                RTopImage = RTopFOV;
                LBotImage = LBotFOV;
                RBotImage = LTopFOV;
                imageFlip = 2;
                imageRot  = -90;
            end
            lengthX = fovY;
            lengthY = fovX;
        elseif(TOrientFOV=='R')
            if(ROrientFOV=='A')
                % rotate counterclock 1
                LTopImage = RTopFOV;
                RTopImage = RBotFOV;
                LBotImage = LTopFOV;
                RBotImage = LBotFOV;
                imageFlip = 0;
                imageRot  = 90;
            else
                % mirror + counterclock 1
                LTopImage  = LTopFOV;
                RTopImage  = LBotFOV;
                LBotImage  = RTopFOV;
                RBotImage  = RBotFOV;
                imageFlip = 2;
                imageRot  = 90;
            end
            lengthX = fovY;
            lengthY = fovX;
        else
            if(ROrientFOV=='R')
                % rotate clockwise 2
                LTopImage = RBotFOV;
                RTopImage = LBotFOV;
                LBotImage = RTopFOV;
                RBotImage = LTopFOV;
                imageFlip = 0;
                imageRot  = 180;
            else
                % mirror + rotate
                LTopImage  = LBotFOV;
                RTopImage  = RBotFOV;
                LBotImage  = LTopFOV;
                RBotImage  = RTopFOV;
                imageFlip = 2;
                imageRot  = 180;
            end
            lengthX = fovX;
            lengthY = fovY;
        end
        
        
    case 'FHLR'
        TOrientImage = 'H';
        BOrientImage = 'F';
        LOrientImage = 'R';
        ROrientImage = 'L';
        
        if(TOrientFOV=='F')
            if(ROrientFOV=='R')
                % rotate clockwise 2
                LTopImage = RBotFOV;
                RTopImage = LBotFOV;
                LBotImage = RTopFOV;
                RBotImage = LTopFOV;
                imageFlip = 0;
                imageRot  = 180;
            else
                % mirror + rotate clockwise 2
                LTopImage = LBotFOV;
                RTopImage = RBotFOV;
                LBotImage = LTopFOV;
                RBotImage = RTopFOV;
                imageFlip = 2;
                imageRot  = 180;
            end
            lengthX = fovX;
            lengthY = fovY;
        elseif(TOrientFOV=='R')
            if(ROrientFOV=='F')
                % mirror + rotate counterclock
                LTopImage  = LBotFOV;
                RTopImage  = RBotFOV;
                LBotImage  = LTopFOV;
                RBotImage  = RTopFOV;
                imageFlip = 2;
                imageRot  = 90;
            else
                % rotate counterclock
                LTopImage = RTopFOV;
                RTopImage = RBotFOV;
                LBotImage = LTopFOV;
                RBotImage = LBotFOV;
                imageFlip = 0;
                imageRot  = 90;
            end
            lengthX = fovY;
            lengthY = fovX;
        elseif(TOrientFOV=='L')
            if(ROrientFOV=='F')
                % rotate clockwise 1
                LTopImage = LBotFOV;
                RTopImage = LTopFOV;
                LBotImage = RBotFOV;
                RBotImage = RTopFOV;
                imageFlip = 0;
                imageRot  = -90;
            else
                % mirror + rotate clockwise 1
                LTopImage  = RBotFOV;
                RTopImage  = RTopFOV;
                LBotImage  = LBotFOV;
                RBotImage  = LTopFOV;
                imageFlip = 2;
                imageRot  = -90;
            end
            lengthX = fovY;
            lengthY = fovX;
        else
            if(ROrientFOV=='R')
                % mirror
                LTopImage = RTopFOV;
                RTopImage = LTopFOV;
                LBotImage = RBotFOV;
                RBotImage = LBotFOV;
                imageFlip = 2;
                imageRot  = 0;
            else
                % proper orientation
                LTopImage = LTopFOV;
                RTopImage = RTopFOV;
                LBotImage = LBotFOV;
                RBotImage = RBotFOV;
                imageFlip = 0;
                imageRot  = 0;
            end
            lengthX = fovX;
            lengthY = fovY;
        end
end

%% get the direction
if contains(plane.foldoverDir,TOrientImage)
    imageFoldoverDir = 'ROW';
else
    imageFoldoverDir = 'COL';
end

%% get extended points
[LTopExtFOV,RTopExtFOV,LBotExtFOV,RBotExtFOV] = ...
    image.plane.findReconPoints( ...
    LTopImage,RTopImage,LBotImage,RBotImage,outerFOVratio);

%% assign to the image
imagePlane.LTop = LTopImage; % Top-Left
imagePlane.RTop = RTopImage; % Top-Right
imagePlane.LBot = LBotImage; % Bottom-Left
imagePlane.RBot = RBotImage; % Bottom-ROrientFOV
% extended FOV points
imagePlane.LTopExtFOV = LTopExtFOV; % Top-Left
imagePlane.RTopExtFOV = RTopExtFOV; % Top-Right
imagePlane.LBotExtFOV = LBotExtFOV; % Bottom-Left
imagePlane.RBotExtFOV = RBotExtFOV; % Bottom-ROrientFOV
% orientation
imagePlane.TOrient = TOrientImage;
imagePlane.BOrient = BOrientImage;
imagePlane.LOrient = LOrientImage;
imagePlane.ROrient = ROrientImage;
% flip and rotation of the image
imagePlane.flip = imageFlip; % 0: no flip / 1: flip first dimension / 2 flip second
imagePlane.rot  = imageRot; % angle in degrees to flip the image (90, -90, 180...)
% dimensions in plane (image FOV)
imagePlane.fovX = fovX;
imagePlane.fovY = fovY;
% foldover direction
imagePlane.foldoverDir = imageFoldoverDir;
% length of image
imagePlane.lengthX = lengthX;
imagePlane.lengthY = lengthY;