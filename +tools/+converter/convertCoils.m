
%% COILS
clear all; close all; clc;

% basePath = 'C:\Users\jorge\CODE\Coil-Maps\CoilData\';
% newPath  = 'C:\Users\jorge\CODE\Coil-Maps\CoilDataNew\';
basePath    = '/efs-mount-point/edt_tool/FILES/coil_models/';
newPath     = '/efs-mount-point/edt_tool/FILES/MODULAR/coil_models/';

coil_list{1}.name = 'OptimalRx';
coil_list{1}.istx = 0;
coil_list{1}.path = [];
coil_list{1}.allowParallelIm        = 0;
coil_list{2}.name = 'GantryBirdcageRx';
coil_list{2}.istx = 1;
coil_list{2}.allowParallelIm        = 0;
coil_list{2}.path = 'Gantry_Birdcage_16Rungs_R350_L600.mat';
coil_list{3}.name = 'HeadBirdcageRx';
coil_list{3}.istx = 1;
coil_list{3}.allowParallelIm        = 0;
coil_list{3}.path = 'Head_Birdcage_16Rungs_R150_L400.mat';
coil_list{4}.name = 'SingleLeftRx';
coil_list{4}.istx = 0;
coil_list{4}.allowParallelIm        = 0;
coil_list{4}.path = 'Head_Array_Wire_R150_1x_Conformal_150x100_Left.mat';
coil_list{5}.name = 'SingleRightRx';
coil_list{5}.istx = 0;
coil_list{5}.allowParallelIm        = 0;
coil_list{5}.path = 'Head_Array_Wire_R150_1x_Conformal_150x100_Right.mat';
coil_list{6}.name = 'SinglePosteriorRx';
coil_list{6}.istx = 0;
coil_list{6}.allowParallelIm        = 0;
coil_list{6}.path = 'Head_Array_Wire_R150_1x_Conformal_150x100_Posterior.mat';
coil_list{7}.name = 'SingleAnteriorRx';
coil_list{7}.istx = 0;
coil_list{7}.allowParallelIm        = 0;
coil_list{7}.path = 'Head_Array_Wire_R150_1x_Conformal_150x100_Anterior.mat';
coil_list{8}.name = '4xConformalRx';
coil_list{8}.istx = 0;
coil_list{8}.allowParallelIm        = 1;
coil_list{8}.parallelImDirection{1} = 'AP';
coil_list{8}.parallelImDirection{2} = 'RL';
coil_list{8}.path = 'Head_Array_Wire_R150_4x_Conformal_150x100.mat';
coil_list{9}.name = '8xConformalRx';
coil_list{9}.istx = 0;
coil_list{9}.allowParallelIm        = 1;
coil_list{9}.parallelImDirection{1} = 'AP';
coil_list{9}.parallelImDirection{2} = 'RL';
coil_list{9}.path = 'Head_Array_Wire_R150_8x_Conformal_150x100.mat';
coil_list{10}.name = '8xPlanarRx';
coil_list{10}.istx = 0;
coil_list{10}.allowParallelIm        = 1;
coil_list{10}.parallelImDirection{1} = 'AP';
coil_list{10}.parallelImDirection{2} = 'RL';
coil_list{10}.path = 'Head_Array_Wire_R150_8x_Planar_150x100.mat';

for ii = 2:length(coil_list)
    
    
    load(sprintf('%s%s',basePath,coil_list{ii}.path));
    
    coilData.freq = coil_data.coil_freq;
    coilData.coilConfig = coil_data.coil_config;
    coilData.numCoils = coil_data.coil_number;
    coilData.coilType = coil_data.coil_type;
    coilData.coilLength = coil_data.coil_length;
    if isfield(coil_data, 'coil_width')
        coilData.coilWidth = coil_data.coil_width;
    else
        coilData.coilWidth = [];
    end
    if isfield(coil_data, 'coil_rungs')
        coilData.coilRungs = coil_data.coil_rungs;
    else
        coilData.coilRungs = [];
    end
    coilData.coilRadius = coil_data.coil_radius;
    coilData.coilRotation = coil_data.coil_rotation;
    coilData.coilIsocenter = coil_data.coil_isocenter;
    coilData.portFeed = coil_data.port_feed;
    coilData.portConfig = coil_data.port_config;
    coilData.coilModel = coil_data.coil_mode;
    coilData.coilTune = coil_data.coil_tune;
    coilData.coilCap = coil_data.coil_cap;
    coilData.geometry = coil_data.coil_geometry;
    coilData.mapsFile = coil_data.maps_file;
    if isfield(coil_data, 'b1m_isocenter')
        coilData.b1mIsocenter = coil_data.b1m_isocenter;
    else
        coilData.b1mIsocenter = 1.0;
    end
    if isfield(coil_data, 'b1p_isocenter')
        coilData.b1pIsocenter = coil_data.b1p_isocenter;
    else
        coilData.b1pIsocenter = 1.0;
    end
    
    
%     coilMaps = load(coilData.mapsFile);
%     maps.bx = coilMaps.bx;
%     maps.by = coilMaps.by;
%     maps.bz = coilMaps.bz;
%     maps.ex = coilMaps.ex;
%     maps.ey = coilMaps.ey;
%     maps.ez = coilMaps.ez;
%     maps.dim = coilMaps.maps_dim;
%     maps.sosB1m = coilMaps.sos_b1m;
%     maps.sosErms = coilMaps.sos_erms;
%     maps.spatial = coilMaps.spatial;
    
    mapsfilename = split(coil_data.maps_file,'/');
    load(sprintf('%s%s',basePath,mapsfilename{end}));
    %bx = bx;
    %by = by;
    %bz = bz;
    %ex = ex;
    %ey = ey;
    %ez = ez;
    dim = maps_dim;
    sosB1m = sos_b1m;
    sosB1p = abs(sos_b1m);
    sosErms = sos_erms;
    %spatial = spatial;

    coilData.mapsFile = sprintf('%s%s',newPath,mapsfilename{end});
    save(sprintf('%s%s',newPath,coil_list{ii}.path), 'coilData', '-v7.3');
    save(coilData.mapsFile, ...
        'bx', 'by', 'bz', 'ex', 'ey', 'ez', ...
        'dim', 'sosB1m', 'sosB1p', 'sosErms', 'spatial', '-v7.3');
    
end
    
