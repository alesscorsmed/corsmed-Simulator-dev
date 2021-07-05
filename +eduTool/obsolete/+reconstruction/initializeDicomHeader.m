function info = initializeDicomHeader(info,dicomInfoNew)

% Update the values
info.InstanceCreationDate    = datestr(now,'yyyyMMdd');
info.InstanceCreationTime    = datestr(now,'HH:MM:SS');
info.SeriesDate              = datestr(now,'yyyyMMdd');
info.AcquisitionDate         = datestr(now,'yyyyMMdd');
info.SeriesTime              = datestr(now,'HH:MM:SS');
info.AcquisitionTime         = datestr(now,'HH:MM:SS');
info.InstitutionName         = 'Corsmed AB';
info.InstitutionAddress      = 'Lund, Sweden';
info.StationName             = 'Educational Tool';
info.StudyDescription        = dicomInfoNew.StudyDescription;  %'Torso Model - Education'
info.SeriesDescription       = dicomInfoNew.SeriesDescription; % 'Gradient Recalled Echo'
info.OperatorsName           = dicomInfoNew.OperatorsName; % 'Christos G. Xanthis'
info.PatientName.FamilyName  = dicomInfoNew.FamilyName; % 'Torso Model'
info.PatientID               = dicomInfoNew.PatientID; % 'v20200316d'
info.PatientBirthDate        = '19901225';
info.PatientSex              = 'M';
info.PatientAge              = '35';
info.PatientSize             = 1.85;
info.PatientWeight           = 90;
info.BodyPartExamined        = dicomInfoNew.BodyPartExamined; % 'Brain'
info.ScanningSequence        = 'GR';
info.NumberOfAverages        = dicomInfoNew.NumberOfAverages;
info.ImagingFrequency        = dicomInfoNew.ImagingFrequency;
info.ImagedNucleus           = '1H';
info.MagneticFieldStrength   = dicomInfoNew.MagneticFieldStrength;
info.NumberOfPhaseEncodingSteps = dicomInfoNew.NumberOfPhaseEncodingSteps;
info.PercentSampling         = 100;
info.PercentPhaseFieldOfView = 100;
info.PixelBandwidth          = dicomInfoNew.PixelBandwidth;
info.SoftwareVersions        = 'v20200316d';
info.ProtocolName            = [];
info.TransmitCoilName        = dicomInfoNew.TransmitCoilName; % 'BODY'
if strcmp(dicomInfoNew.foldoverDir,'ROW')
    info.AcquisitionMatrix	= [0;dicomInfoNew.NumberOfFreqEncodingSteps;dicomInfoNew.NumberOfPhaseEncodingSteps;0];
    info.PixelSpacing       = flipud(dicomInfoNew.PixelSpacing);
    info.Height             = dicomInfoNew.NumberOfPhaseEncodingSteps;
    info.Width              = dicomInfoNew.NumberOfFreqEncodingSteps;
else
    info.AcquisitionMatrix	= [dicomInfoNew.NumberOfFreqEncodingSteps;0;0;dicomInfoNew.NumberOfPhaseEncodingSteps];
    info.PixelSpacing       = dicomInfoNew.PixelSpacing;
    info.Height             = dicomInfoNew.NumberOfFreqEncodingSteps;
    info.Width              = dicomInfoNew.NumberOfPhaseEncodingSteps;
end
info.InPlanePhaseEncodingDirection = dicomInfoNew.foldoverDir; % 'ROW' or 'COL'
info.FlipAngle               = dicomInfoNew.FlipAngle;
info.VariableFlipAngleFlag   = 'N'; % 'Y' or 'N'
info.SAR                     = [];
info.PatientPosition         = 'HFS';
info.AcquisitionNumber       = dicomInfoNew.AcquisitionNumber;
info.ImagePositionPatient    = dicomInfoNew.ImagePositionPatient; % [0;192;0]
info.ImageOrientationPatient = dicomInfoNew.ImageOrientationPatient; % [1;0;0;0;1;0]
info.SliceLocation           = dicomInfoNew.SliceLocation;
info.ImageComments           = dicomInfoNew.ImageComments; % 'My comments'
info.StudyID                 = dicomInfoNew.StudyID; % 'StudyID: 1'
