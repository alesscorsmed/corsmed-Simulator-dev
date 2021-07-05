function generateDicomFile(sliceContrast, slicePlane, ...
    bodyPartName, dicomInfo, dicomStructInitial)
%
% IMAGE.DICOM.generateDicomFile
%
%     generates the dicom file
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================

%% Update info to be used in the dicom header
dicomInfo.BodyPartExamined          	= ['Corsmed - ',bodyPartName];
dicomInfo.AcquisitionNumber            	= sliceContrast.uniqueID;
% Find the orientation and position of the imaging plane
[dicomInfo] = image.dicom.findDicomOrientationPosition( ...
    slicePlane, dicomInfo);

%% Formulate the dicom header
dicomstruct = image.dicom.initializeDicomHeader( ...
    dicomStructInitial, dicomInfo);

%% Generate the .dcm file
dicomImageMatrix = abs(sliceContrast.image);
dicomwrite(dicomImageMatrix, sliceContrast.dcmName, dicomstruct,...
    'CreateMode', 'copy');
disp('.dcm file created!')
