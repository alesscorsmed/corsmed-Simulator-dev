function dicfile = createDicom(inputImage, dicomstruct, exportDicomName, savePath)
%The function below takes as input:
%>>> The Image we want to form as dicom in 2D array format or in 
%file format (" or ' quotes)
%>>> And the struct that we want our dicom image to have as a Header
%>>> The name we want to give to the exported dicom file in (''or "") format
%>>> Optionally give the path we want the dicom file to be saved in. (''or
%"") format. If we dont want a specific path just add empty char ''.
%
% dicomstruct holds the dicom tags for the specific dicom file. Initially,
% the dicomstruct is loaded from the emptyDicomStruct.mat file using the
% following command:
% dicomstruct = load('emptyDicomStruct.mat');
% The emptyDicomStruct.mat was created using the following command where
% aaa is the struct that holds all the dicom tags:
% save('emptyDicomStruct.mat','-struct','aaa')

if(isstring(inputImage))
    inputImage = convertStringsToChars(inputImage);
end

if(ischar(inputImage))
    inputImage = imread(inputImage);    
end

if(isempty(savePath) || strlength(savePath) <= 1)
    savePath = pwd;
end

dicfile = fullfile(savePath,sprintf('%s.dcm',exportDicomName));

if(length(dicomstruct.Filename) <=1 ) 
    dicomwrite(inputImage,dicfile);
else
    dicomwrite(inputImage,dicfile,dicomstruct);
end