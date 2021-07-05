function [imageData] = defineContrastProcess(imageData, seqFamilyName)
%
% IMAGE.DICOM.DEFINECONTRASTPROCESSIMAGE
%
% Defines the actual number of contrasts of the image 
% and how to process the recon images to generate them
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

switch lower(seqFamilyName)
    case 'gre-oop'
        % define how to process the contrasts
        if imageData.b0map
            imageData.numContrasts  = imageData.numContrasts + 3;
        else
            imageData.numContrasts  = 2*imageData.numContrasts;
        end
        imageData.processContrast   = 'waterFatOOP';
    case 'swi-2d-gre'
        % define how to process the contrasts
        imageData.numContrasts      = imageData.numContrasts + 2;
        imageData.processContrast   = 'swi';
    case 'molli'
        % define how to process the contrasts
        imageData.numContrasts      = imageData.numContrasts + 1;
        imageData.processContrast   = 'molli';  
    otherwise
        imageData.processContrast   = 'none';
end