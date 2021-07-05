function info = createDefaultDicomHeader()

% Default dicom header created by MATLAB
dicomwrite(mat2gray(0.5*ones(162,192)),'test.dcm',...
    'ObjectType','MR Image Storage');
infoDefault     = dicominfo('test.dcm');

% Add the following fields
infoDefault.InstanceCreationDate    = [];
infoDefault.InstanceCreationTime    = [];
infoDefault.SeriesDate              = [];
infoDefault.AcquisitionDate         = [];
infoDefault.SeriesTime              = [];
infoDefault.AcquisitionTime         = [];
infoDefault.InstitutionName         = [];
infoDefault.InstitutionAddress      = [];
infoDefault.StationName             = [];
infoDefault.StudyDescription        = [];
infoDefault.SeriesDescription       = [];
infoDefault.OperatorsName           = [];
infoDefault.PatientAge              = [];
infoDefault.PatientSize             = [];
infoDefault.PatientWeight           = [];
infoDefault.BodyPartExamined        = [];
infoDefault.ScanningSequence        = [];
infoDefault.NumberOfAverages        = [];
infoDefault.ImagingFrequency        = [];
infoDefault.ImagedNucleus           = []; 
infoDefault.MagneticFieldStrength   = [];
infoDefault.NumberOfPhaseEncodingSteps = [];
infoDefault.PercentSampling         = [];
infoDefault.PercentPhaseFieldOfView = [];
infoDefault.PixelBandwidth          = [];
infoDefault.SoftwareVersions        = [];
infoDefault.ProtocolName            = [];
infoDefault.TransmitCoilName        = [];
infoDefault.AcquisitionMatrix       = [];
infoDefault.InPlanePhaseEncodingDirection = [];
infoDefault.FlipAngle               = [];
infoDefault.VariableFlipAngleFlag   = [];
infoDefault.SAR                     = [];
infoDefault.PatientPosition         = [];
infoDefault.AcquisitionNumber       = [];
infoDefault.ImagePositionPatient    = [];
infoDefault.ImageOrientationPatient = [];
infoDefault.SliceLocation           = [];
infoDefault.ImageComments           = [];
infoDefault.PixelSpacing            = [];

% Dicom header from a .dcm file generated by the MR scanner
info            = dicominfo(['inputs',filesep,'1605.dcm']);
trueDicomfields = fieldnames(info);

for i=1:size(trueDicomfields,1)
    if ~isfield(infoDefault,trueDicomfields{i})
        info = rmfield(info,trueDicomfields{i});
    else
        evalc(['info.',trueDicomfields{i},'=infoDefault.',trueDicomfields{i},';']);
    end
end