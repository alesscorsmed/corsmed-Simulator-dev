function [dicomInfo] = generateDicomInfo(acquisition, mrSystem, expControl)
%
% IMAGE.DICOM.GENERATEDICOMINFO
%
%     generates dicom info for file
%
% INPUT
%   acquisition        structure with acquisition data
%   mrSystem           mrSystem info
%   expControl         experiment control data, with commLocalDB
%
% OUTPUT
%   dicomInfo          initialized struct
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'image.dicom.generateDicomInfo';
if (nargin < 3)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

% Update the dicomInfo structure
dicomInfo.StudyDescription              = 'Corsmed - Educational Tool';
dicomInfo.SeriesDescription             = acquisition.data.pulseSeqFamilyName;
dicomInfo.OperatorsName                 = ['User-',num2str(expControl.userID)];
dicomInfo.FamilyName                    = 'Corsmed';
dicomInfo.PatientID                     = expControl.versionNum;
dicomInfo.BodyPartExamined          	= []; % TO BE FILLED
dicomInfo.NumberOfAverages            	= acquisition.data.NEX;
dicomInfo.ImagingFrequency            	= mrSystem.b0*42.56;
dicomInfo.MagneticFieldStrength        	= mrSystem.b0;
dicomInfo.NumberOfPhaseEncodingSteps   	= acquisition.data.numPE;
dicomInfo.NumberOfFreqEncodingSteps    	= acquisition.data.numFE;
dicomInfo.PixelBandwidth               	= acquisition.data.rxBW/...
    acquisition.data.numFE;
dicomInfo.TransmitCoilName             	= 'BODY';
dicomInfo.FlipAngle                   	= acquisition.mainRF.flipAngle;
dicomInfo.AcquisitionNumber            	= 0; % TO BE FILLED
dicomInfo.SliceLocation               	= 0; % Usually this attribute is not needed
dicomInfo.StudyID                      	= num2str(acquisition.data.pulseSeqNum);
dicomInfo.ImageComments               	= 'None';   % FIX IT
dicomInfo.foldoverDir                  	= ''; % TO BE FILLED
if contains(lower(acquisition.data.pulseSeqFamilyName), 'gre-oop')
    dicomInfo.EchoTime                       = 2*acquisition.data.deltaTimeOOP;
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done ', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
