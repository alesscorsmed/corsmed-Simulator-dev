function sessionData = findAnatomicalModelPath(sessionData,courseId)
%
% EDUTOOL.CONNECTDATABASE
%
%	returns information for the anatomical Model that will be used
%
% INPUT
%
%   applicationStruct         struct with information for the start of an exper.
%   courseId                 the unique id representing the course
%
% OUTPUT
%
%   applicationStruct        struct with Data necessary for start of Edutool
%
%========================  CORSMED AB © 2020 ==============================
%

%Choose Model
switch str2double(courseId)
    
    case 4        
        fprintf('XCAT anatomical model ')
        modelName = 'Torso Model';
        defaultAnatModelPath = '/efs-mount-point/edt_tool/FILES/anatomical_models/anatomical_model_20190326_234501.mat';
        defaultCoilsModelPath   = '';
        PDfluctuationFactor     = 0.2;
        default_gridStep        = [0.001,0.001,0.001];
        fatTissuesIDs           = [];
        PDinhomPathFile         = '/efs-mount-point/edt_tool/FILES/anatomical_models/20200506_XCAT_PDinhom.mat';
        volumeGridSize          = [];
    case 6     
        fprintf('MIDA anatomical model ')
        modelName = 'MIDA';
    %     defaultAnatModelPath = '/efs-mount-point/edt_tool/FILES/anatomical_models/MIDA_coremri_20191201.mat';
    %     defaultAnatModelPath = '/efs-mount-point/edt_tool/FILES/anatomical_models/MIDA_coremri_20200202.mat';
    %     defaultAnatModelPath = '/efs-mount-point/edt_tool/FILES/anatomical_models/MIDA_coremri_20200210_coilmaps.mat';
    %     defaultAnatModelPath = '/efs-mount-point/edt_tool/FILES/anatomical_models/MIDA_coremri_20200212_birdcageANDsingleCoil.mat';
        defaultAnatModelPath    = '/efs-mount-point/edt_tool/FILES/anatomical_models/MIDA_coremri_20200322_AnatomicalModel.mat';
        defaultCoilsModelPath   = '/efs-mount-point/edt_tool/FILES/anatomical_models/MIDA_coremri_20200322_Coils.mat';
        PDfluctuationFactor     = 0.05;
        default_gridStep        = [0.0005,0.0005,0.0005];
        fatTissuesIDs           = [43,62];
        PDinhomPathFile         = '/efs-mount-point/edt_tool/FILES/anatomical_models/20200506_MIDA_PDinhom.mat';
        volumeGridSize          = [480,350,480];
    case 7
        fprintf('Five cylinders model  ')
        modelName               = 'Cylinders';
        defaultAnatModelPath    = '/efs-mount-point/edt_tool/FILES/anatomical_models/20200320_fiveCylinders_TwoInnerCylinders.mat';
        defaultCoilsModelPath   = '';
        PDfluctuationFactor     = 0.05;
        default_gridStep        = [0.001,0.001,0.001];
        fatTissuesIDs           = [];
        PDinhomPathFile         = '';
        volumeGridSize          = [];
    case 8
        fprintf('Voxelman): ')
        modelName               = 'Voxelman';
        defaultAnatModelPath    = '/efs-mount-point/edt_tool/FILES/anatomical_models/20200520d_VoxelManModel4EduTool_withAir.mat';
        defaultCoilsModelPath   = '';
        PDfluctuationFactor     = 0.05;
        default_gridStep        = [0.001,0.001,0.001];
        fatTissuesIDs           = [];
        PDinhomPathFile         = '';
        volumeGridSize          = [];
end

sessionData.anatModel.modelName               = modelName;
sessionData.anatModel.defaultAnatModelPath    = defaultAnatModelPath;
sessionData.anatModel.defaultCoilsModelPath   = defaultCoilsModelPath;
sessionData.anatModel.PDfluctuationFactor     = PDfluctuationFactor;
sessionData.anatModel.default_gridStep        = default_gridStep;
sessionData.anatModel.fatTissuesIDs           = fatTissuesIDs;
sessionData.anatModel.PDinhomPathFile         = PDinhomPathFile;
sessionData.anatModel.volumeGridSize          = volumeGridSize;
