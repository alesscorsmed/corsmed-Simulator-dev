clear all; close all; clc;

experimentList = dir('/efs-mount-point-MATLAB/EduTool-Jorge/TESTV2/V2*.json');

for ii=1:length(experimentList)
    
    % get the names of the experiments
    testName        = experimentList(ii).name;
    v2FileName      = sprintf('%s/%s',experimentList(ii).folder,testName);
    
    testNameV1      = testName;
    testNameV1(2)   = '1';
    v1FileName      = sprintf('%s/%s',experimentList(ii).folder,testNameV1);
    
    errorName       = ['Error', testName(3:end)];
    errorFileName      = sprintf('%s/%s',experimentList(ii).folder,errorName);
    
    
    experimentDataFile = sprintf(v1FileName);
    fid = fopen(experimentDataFile,'r');
    experimentDataV1 = jsondecode(fread(fid,inf,'*char').');
    fclose(fid);
    
    experimentDataFile = sprintf(v2FileName);
    fid = fopen(experimentDataFile,'r');
    experimentDataV2 = jsondecode(fread(fid,inf,'*char').');
    fclose(fid);
    
    errorStruct = [];
    
    %% expControl
    
    % model
    errorStruct.expControl.model = [];
    [errorStruct.expControl.model] = tools.misc.deepCompareStruct(...
        experimentDataV2.expControl.model,...
        experimentDataV1.expControl.model,...
        errorStruct.expControl.model);
    
    % simulation
    errorStruct.expControl.simulation = [];
    [errorStruct.expControl.simulation] = tools.misc.deepCompareStruct(...
        experimentDataV2.expControl.simulation,...
        experimentDataV1.expControl.simulation,...
        errorStruct.expControl.simulation);
    
    % sequence
    errorStruct.expControl.sequence = [];
    [errorStruct.expControl.sequence] = tools.misc.deepCompareStruct(...
        experimentDataV2.expControl.sequence,...
        experimentDataV1.expControl.sequence,...
        errorStruct.expControl.sequence);
    
    % mrsystem
    errorStruct.expControl.mrSystem = [];
    [errorStruct.expControl.mrSystem] = tools.misc.deepCompareStruct(...
        experimentDataV2.expControl.mrSystem,...
        experimentDataV1.expControl.mrSystem,...
        errorStruct.expControl.mrSystem);
    
    % motionSpecs
    errorStruct.expControl.motionSpecs = [];
    [errorStruct.expControl.motionSpecs] = tools.misc.deepCompareStruct(...
        experimentDataV2.expControl.motionSpecs,...
        experimentDataV1.expControl.motionSpecs,...
        errorStruct.expControl.motionSpecs);
    
    %% acquisition
    
    % data
    errorStruct.acquisition.data = [];
    [errorStruct.acquisition.data] = tools.misc.deepCompareStruct(...
        experimentDataV2.acquisition.data,...
        experimentDataV1.acquisition.data,...
        errorStruct.acquisition.data);
    
    % mainRF
    errorStruct.acquisition.mainRF = [];
    [errorStruct.acquisition.mainRF] = tools.misc.deepCompareStruct(...
        experimentDataV2.acquisition.mainRF,...
        experimentDataV1.acquisition.mainRF,...
        errorStruct.acquisition.mainRF);
    
    % refRF
    errorStruct.acquisition.refRF = [];
    [errorStruct.acquisition.refRF] = tools.misc.deepCompareStruct(...
        experimentDataV2.acquisition.refRF,...
        experimentDataV1.acquisition.refRF,...
        errorStruct.acquisition.refRF);
    
    % prepIR
    errorStruct.acquisition.prepIR = [];
    [errorStruct.acquisition.prepIR] = tools.misc.deepCompareStruct(...
        experimentDataV2.acquisition.prepIR,...
        experimentDataV1.acquisition.prepIR,...
        errorStruct.acquisition.prepIR);
    
    % encPG
    errorStruct.acquisition.encPG = [];
    [errorStruct.acquisition.encPG] = tools.misc.deepCompareStruct(...
        experimentDataV2.acquisition.encPG,...
        experimentDataV1.acquisition.encPG,...
        errorStruct.acquisition.encPG);
    
    % noise
    errorStruct.acquisition.noise = [];
    [errorStruct.acquisition.noise] = tools.misc.deepCompareStruct(...
        experimentDataV2.acquisition.noise,...
        experimentDataV1.acquisition.noise,...
        errorStruct.acquisition.noise);
    
    errorStructData = jsonencode(errorStruct);
    errorDataFile = sprintf(errorFileName);
    fid = fopen(errorDataFile,'w');
    fwrite(fid, errorStructData, 'char');
    fclose(fid);
    
end
