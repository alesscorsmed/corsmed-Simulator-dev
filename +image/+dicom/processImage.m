function [contrastData] = processImage(contrastData, ...
imageMap, tissueMask, slicePlane)
%
% IMAGE.DICOM.PROCESSIMAGE
%
% The normalizedImageResized is the resized image so as to keep the correct
% ratio between the vertical and horizontal axis of the image based on the 
% size of FOV and kspace. The normalizedImage is the original image where
% no "stretch" has been applied.
%
% Flipping and rotation is also applied
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

%functionName = 'image.dicom.processImage';

% Image size
nX = slicePlane.sizeX;
nY = slicePlane.sizeY;
% get FOV of image
fovX = slicePlane.fovX;
fovY = slicePlane.fovY;
% pixel size
pixelSizeX = fovX/nX;
pixelSizeY = fovY/nY;

% resize
if fovX<fovY
    resizeFactor = (nX/nY)*fovY/fovX;
    resizeFactorArray = [resizeFactor*nY nX];
elseif fovX>fovY
    resizeFactor = (nY/nX)*fovX/fovY;
    resizeFactorArray = [nY resizeFactor*nX];
else
    if pixelSizeY>pixelSizeX
        resizeFactor = (nX/nY)*fovY/fovX;
        resizeFactorArray = [resizeFactor*nY nX];
    elseif pixelSizeY<pixelSizeX
        resizeFactor = (nY/nX)*fovX/fovY;
        resizeFactorArray = [nY resizeFactor*nX];
    else
        resizeFactorArray = [nY nX];
    end
end

%% normalize image
contrastImage   = imageMap.';
contrastMask    = tissueMask.';

%% resize image and mask
contrastImage = imresize(contrastImage ,resizeFactorArray);
contrastMask  = imresize(contrastMask  ,resizeFactorArray, 'nearest');

%% correct orientation
if (slicePlane.flip > 0)
    contrastImage = flip(contrastImage,slicePlane.flip);
    contrastMask  = flip(contrastMask ,slicePlane.flip);
end
contrastImage = imrotate(contrastImage,slicePlane.rot);
contrastMask  = imrotate(contrastMask ,slicePlane.rot);

%% assign
contrastData.normalizedImage    = mat2gray(abs(contrastImage));
contrastData.resizeFactorArray  = resizeFactorArray;
contrastData.image              = contrastImage;
contrastData.mask               = contrastMask;
