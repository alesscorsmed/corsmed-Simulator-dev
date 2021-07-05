function generateDicomFile(outputImage,acquisition,spinModel,...
    expControl,sessionData,mrSystem,dicomStructInitial)

%% Update info to be used in the dicom header
% Find the orientation and position of the imaging plane
[dicomInfo.ImagePositionPatient,dicomInfo.ImageOrientationPatient,...
    dicomInfo.PixelSpacing] = ...
    eduTool.reconstruction.findDicomOrientationPosition(outputImage);

% Update the dicomInfo structure
dicomInfo.StudyDescription              = 'Corsmed - Educational Tool';
dicomInfo.SeriesDescription             = acquisition.data.pulseSeqFamilyName;
dicomInfo.OperatorsName                 = ['User-',num2str(expControl.userID)];
dicomInfo.FamilyName                    = 'Corsmed';
dicomInfo.PatientID                     = sessionData.versionNum;
dicomInfo.BodyPartExamined          	= ['Corsmed - ',spinModel.bodyPartName];
dicomInfo.NumberOfAverages            	= acquisition.data.NEX;
dicomInfo.ImagingFrequency            	= mrSystem.b0*42.56;
dicomInfo.MagneticFieldStrength        	= mrSystem.b0;
dicomInfo.NumberOfPhaseEncodingSteps   	= acquisition.data.numPE;
dicomInfo.NumberOfFreqEncodingSteps    	= acquisition.data.numFE;
dicomInfo.PixelBandwidth               	= acquisition.data.rxBW/...
    acquisition.data.numFE;
dicomInfo.TransmitCoilName             	= 'BODY';
dicomInfo.FlipAngle                   	= acquisition.mainRF.flipAngle;
dicomInfo.AcquisitionNumber            	= outputImage.uniqueID;
dicomInfo.SliceLocation               	= 0; % Usually this attribute is not needed
dicomInfo.StudyID                      	= num2str(acquisition.data.pulseSeqNum);
dicomInfo.ImageComments               	= 'None';   % FIX IT
dicomInfo.foldoverDir                  	= outputImage.foldoverDir;

%% Formulate the dicom header
dicomstruct = eduTool.reconstruction.initializeDicomHeader(dicomStructInitial,...
    dicomInfo);

%% Generate the .dcm file
dicomwrite(outputImage.image,outputImage.dcmName,dicomstruct,...
    'CreateMode','copy')
disp('.dcm file created!')