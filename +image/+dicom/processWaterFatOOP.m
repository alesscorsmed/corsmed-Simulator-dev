function [contrastSpace,contrastText,isContrast] = processWaterFatOOP(...
    multiContrastSlice, imageData)
%
% IMAGE.DICOM.PROCESSWATERFATOOP
%
% Process the 2-echo contrasts to generate the 4 contrasts:
%   Echo 1 image
%   Echo 2 image
%   Water  image
%   Fat    image
%
% If B0map == yes (through the UI), an extra echo is acquired and the B0map
% is generated
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'image.dicom.processWaterFatOOP';

%% allocate
[nX,nY,nC] = size(multiContrastSlice);
% assert(nC == 2, ...
%     sprintf('error: %s number of contrast %d different from 2\n',...
%     functionName,nC));

if (imageData.b0map == 1)
   deltaTE = imageData.dicomInfo.EchoTime/2;
   contrastSpace = zeros(nX,nY,6);
   contrastSpace(:,:,1) = multiContrastSlice(:,:,1);
   contrastSpace(:,:,2) = multiContrastSlice(:,:,2);
   contrastSpace(:,:,3) = multiContrastSlice(:,:,3);
   contrastSpace(:,:,4) = abs((multiContrastSlice(:,:,1) + multiContrastSlice(:,:,2))/2);
   contrastSpace(:,:,5) = abs((multiContrastSlice(:,:,1) - multiContrastSlice(:,:,2))/2);
   contrastSpace(:,:,6) = angle(multiContrastSlice(:,:,3) ./ multiContrastSlice(:,:,1))/(2*pi*2*deltaTE);
   
   %% text
   contrastText{1} = 'Echo 1';
   contrastText{2} = 'Echo 2';
   contrastText{3} = 'Echo 3';
   contrastText{4} = 'Water (OOP)';
   contrastText{5} = 'Fat   (OOP)';
   contrastText{6} = 'Off-res map';
   
   %% is contrast or map
   isContrast       = ones(6,1);
   isContrast(6)    = 0;
   
else
    contrastSpace = zeros(nX,nY,4);
    %% compute the corresponding images
    contrastSpace(:,:,1) = multiContrastSlice(:,:,1);
    contrastSpace(:,:,2) = multiContrastSlice(:,:,2);
    %% for the IN and OUT images are done with the magnitude
    %  in general cases, complex need to be taken into account
    %  in our (ideal) case, contrast are generated from magnitude
    contrastSpace(:,:,3) = abs((multiContrastSlice(:,:,1) + multiContrastSlice(:,:,2))/2);
    contrastSpace(:,:,4) = abs((multiContrastSlice(:,:,1) - multiContrastSlice(:,:,2))/2);
    
    %% text
    contrastText{1} = 'Echo 1';
    contrastText{2} = 'Echo 2';
    contrastText{3} = 'Water (OOP)';
    contrastText{4} = 'Fat   (OOP)';
    
    %% is contrast or map
    isContrast       = ones(4,1);
end
