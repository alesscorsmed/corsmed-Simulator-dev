function [outputImages] = resizeImage(iMap,reconstruction)
% The normalizedImageResized is the resized image so as to keep the correct
% ratio between the vertical and horizontal axis of the image based on the 
% size of FOV and kspace. The normalizedImage is the original image where
% no "stretch" has been applied.

FOV(1,1)    = reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.x;
FOV(1,2)    = reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.y;

reconNx     = reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.x;
reconNy     = reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.y;

outputImages.normalizedImage = mat2gray(abs(iMap));

pixelx = FOV(1,1)/reconNx;
pixely = FOV(1,2)/reconNy;

if FOV(1,1)<FOV(1,2)
    resizeFactor = (reconNx/reconNy)*FOV(1,2)/FOV(1,1);
    outputImages.resizeFactorArray = [resizeFactor*reconNy reconNx];
    outputImages.normalizedImageResized = ...
        imresize(outputImages.normalizedImage,outputImages.resizeFactorArray);
elseif FOV(1,1)>FOV(1,2)
    resizeFactor = (reconNy/reconNx)*FOV(1,1)/FOV(1,2);
    outputImages.resizeFactorArray = [reconNy resizeFactor*reconNx];
    outputImages.normalizedImageResized = ...
        imresize(outputImages.normalizedImage,outputImages.resizeFactorArray);
else
    if pixely>pixelx
        resizeFactor = (reconNx/reconNy)*FOV(1,2)/FOV(1,1);
        outputImages.resizeFactorArray = [resizeFactor*reconNy reconNx];
        outputImages.normalizedImageResized = ...
            imresize(outputImages.normalizedImage,outputImages.resizeFactorArray);
    elseif pixely<pixelx
        resizeFactor = (reconNy/reconNx)*FOV(1,1)/FOV(1,2);
        outputImages.resizeFactorArray = [reconNy resizeFactor*reconNx];
        outputImages.normalizedImageResized = ...
            imresize(outputImages.normalizedImage,outputImages.resizeFactorArray);
    else
        outputImages.resizeFactorArray = [reconNy reconNx];
        outputImages.normalizedImageResized = outputImages.normalizedImage;
    end
end