function [coilSystem] = initializeCoils(courseID, basePath)
%
% DATA.COILSYSTEM.INITIALIZE
%
%	Function that initializes the coilSystem data structure.
%   Returns a coilSystem with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   courseID     ID number of the course running
%
% OUTPUT
%   coilSystem   structure with coil models and info
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'data.coilModel.initialize';
if (nargin < 1)
    courseID = 6; % BrainApp as default
end
if (nargin < 2) || isempty(basePath)
    ME = MException('Model:wrongPath',...
        '%s : Coil model path not available',functionName);
    throw(ME);
end

%% report start
fid = 1;
tTotal = tic();
fprintf(fid, '\n%s : start', functionName);
fprintf(fid, '\n  Loading data from %s', basePath);

%% create list with coils available for the course
switch courseID
    
    case 4
        %% XCAT model
        % default isocenter and transmitter
        systemIsocenter     = [0.5280, 0.5035, 1.4000];
        defaultCoil         = 'GantryBirdcageRx';
        % optimal coil (ideal, no maps)
        coilList{1}.name    = 'OptimalRx';
        coilList{1}.isTx    = 0;
        coilList{1}.parAP   = 0;
        coilList{1}.parRL   = 0;
        coilList{1}.parFH   = 0;
        coilList{1}.file    = [];
        % gantry
        coilList{2}.name    = 'GantryBirdcageRx';
        coilList{2}.isTx    = 1;
        coilList{2}.parAP   = 0;
        coilList{2}.parRL   = 0;
        coilList{2}.parFH   = 0;
        coilList{2}.file    = 'Gantry_Birdcage_16Rungs_R350_L600.mat';

    case 6 
        %% MIDA
        % default isocenter and transmitter
        systemIsocenter     =  [0.0875, 0.1200, 0.1200];
        defaultCoil         = 'GantryBirdcageRx';
        % optimal coil (ideal, no maps)
        coilList{1}.name    = 'OptimalRx';
        coilList{1}.isTx    = 0;
        coilList{1}.parAP   = 0;
        coilList{1}.parRL   = 0;
        coilList{1}.parFH   = 0;
        coilList{1}.file    = [];
        % gantry
        coilList{2}.name    = 'GantryBirdcageRx';
        coilList{2}.isTx    = 1;
        coilList{2}.parAP   = 0;
        coilList{2}.parRL   = 0;
        coilList{2}.parFH   = 0;
        coilList{2}.file    = 'Gantry_Birdcage_16Rungs_R350_L600.mat';
        % head birdcage
        coilList{3}.name    = 'HeadBirdcageRx';
        coilList{3}.isTx    = 1;
        coilList{3}.parAP   = 0;
        coilList{3}.parRL   = 0;
        coilList{3}.parFH   = 0;
        coilList{3}.file    = 'Head_Birdcage_16Rungs_R150_L400.mat';
        % left coil
        coilList{4}.name    = 'SingleLeftRx';
        coilList{4}.isTx    = 0;
        coilList{4}.parAP   = 0;
        coilList{4}.parRL   = 0;
        coilList{4}.parFH   = 0;
        coilList{4}.file    = 'Head_Array_Wire_R150_1x_Conformal_150x100_Left.mat';
        % head birdcage
        coilList{5}.name    = 'SingleRightRx';
        coilList{5}.isTx    = 0;
        coilList{5}.parAP   = 0;
        coilList{5}.parRL   = 0;
        coilList{5}.parFH   = 0;
        coilList{5}.file    = 'Head_Array_Wire_R150_1x_Conformal_150x100_Right.mat';
        % head birdcage
        coilList{6}.name    = 'SinglePosteriorRx';
        coilList{6}.isTx    = 0;
        coilList{6}.parAP   = 0;
        coilList{6}.parRL   = 0;
        coilList{6}.parFH   = 0;
        coilList{6}.file    = 'Head_Array_Wire_R150_1x_Conformal_150x100_Posterior.mat';
        % head birdcage
        coilList{7}.name    = 'SingleAnteriorRx';
        coilList{7}.isTx    = 0;
        coilList{7}.parAP   = 0;
        coilList{7}.parRL   = 0;
        coilList{7}.parFH   = 0;
        coilList{7}.file    = 'Head_Array_Wire_R150_1x_Conformal_150x100_Anterior.mat';
        % head birdcage
        coilList{8}.name    = '4xConformalRx';
        coilList{8}.isTx    = 0;
        coilList{8}.parAP   = 4;
        coilList{8}.parRL   = 4;
        coilList{8}.parFH   = 0;
        coilList{8}.file    = 'Head_Array_Wire_R150_4x_Conformal_150x100.mat';
        % head birdcage
        coilList{9}.name    = '8xConformalRx';
        coilList{9}.isTx    = 0;
        coilList{9}.parAP   = 4;
        coilList{9}.parRL   = 4;
        coilList{9}.parFH   = 0;
        coilList{9}.file    = 'Head_Array_Wire_R150_8x_Conformal_150x100.mat';
        % head birdcage
        coilList{10}.name   = '8xPlanarRx';
        coilList{10}.isTx   = 0;
        coilList{10}.parAP  = 4;
        coilList{10}.parRL  = 4;
        coilList{10}.parFH  = 0;
        coilList{10}.file   = 'Head_Array_Wire_R150_8x_Planar_150x100.mat';

    otherwise
        %% Other cases
        % default isocenter and transmitter
        systemIsocenter     =  [0.0, 0.0, 0.0];
        defaultCoil         = 'OptimalRx';
        % optimal coil (ideal, no maps)
        coilList{1}.name    = 'OptimalRx';
        coilList{1}.isTx    = 1;
        coilList{1}.parAP   = 0;
        coilList{1}.parRL   = 0;
        coilList{1}.parFH   = 0;
        coilList{1}.file    = [];
end
    

%% loop on the coil list and load the data
numModels = 0;
defaultIndex = 1; % by default optimal: will change if defined

for cc = 1:length(coilList)

    numModels = numModels+1;
    try
        load(sprintf('%s/%s', basePath, coilList{cc}.file), 'coilData');
        % % in case we want to change the path
        % [~,mapFile,mapExt] = fileparts(coilData.mapsFile);
        % coilData.mapsFile = sprintf('%s%s',mapFile,mapExt);
        % save(sprintf('%s/%s', basePath, coilList{cc}.file), 'coilData', '-v7.3');
        % assign the coil data
        data = coilData;
        % update the path to the field maps
        data.mapsFile = sprintf('%s/%s',basePath,data.mapsFile);
    catch
        % coil does not have data, treat it as optimal
        data.numCoils       = 1;
        data.coilConfig     = 'optimal';
        data.b1mIsocenter   = 1.0;
        data.b1pIsocenter   = 1.0;
    end
    % general data
    data.name           = coilList{cc}.name; % name
    data.coilIsocenter  = [0.0, 0.0, 0.0]; % isocenter of coil
    data.isTx           = coilList{cc}.isTx; % if coil is also transmit (transceiver coil)
    data.parallelAP     = coilList{cc}.parAP; % max parallelism in AP direction (0 is no parallel)
    data.parallelRL     = coilList{cc}.parRL; % max parallelism in RL direction (0 is no parallel)
    data.parallelFH     = coilList{cc}.parFH; % max parallelism in FH direction (0 is no parallel)
    % default transmitter
    if strcmpi(defaultCoil, data.name)
        defaultIndex = cc;
        b1mScaling     = data.b1mIsocenter;
        b1pScaling     = data.b1pIsocenter;
    end

    if strcmpi(data.coilConfig,'optimal')
        % optimal coil, empty maps
        maps.dim        = [];  % dimensions nx, ny, nz
        maps.spatial    = [];  % spatial coordinates
        maps.b1mSOS     = [];  % b1m Sum-of-Squares map
        maps.b1pSOS     = [];  % b1p Sum-of-Squares map
        maps.bx         = [];  % bx component maps
        maps.by         = [];  % by component maps
        
        % sar data and map
        sar.estSAR = 'N/A';
        sar.avgSAR = 'N/A';
        sar.cm3SAR = 'N/A';
        sar.g10SAR = 'N/A';
        sar.voxSAR = 'N/A';
        sar.mapSAR = [];
        
    else
        % maps data: allocate, do not load yet (load on demmand)
        % maps.dim        = [];  % dimensions nx, ny, nz
        % maps.spatial    = [];  % spatial coordinates
        % maps.b1mSOS     = [];  % b1m Sum-of-Squares map
        % maps.b1pSOS     = [];  % b1p Sum-of-Squares map
        % maps.bx         = [];  % bx component maps
        % maps.by         = [];  % by component maps
        maps = [];
        
        % sar data and map
        % sar.estSAR = [];
        % sar.avgSAR = [];
        % sar.cm3SAR = [];
        % sar.g10SAR = [];
        % sar.voxSAR = [];
        % sar.mapSAR = [];
        sar = [];
    end

    % assign to struct
    coilSystem.coilModel{cc}.data = data;
    coilSystem.coilModel{cc}.maps = maps;
    coilSystem.coilModel{cc}.sar  = sar; 
    
    fprintf(fid, '\n  %20s coil loaded -- parallelism AP=%d/RL=%d/FH=%d',...
        data.name, data.parallelAP, data.parallelRL, data.parallelFH);
    
end

%% in case there is no list, force optimal coil
if (numModels == 0)
    
    numModels = numModels+1;
    
    % general data
    data.name           = 'OptimalRx'; % name
    data.coilIsocenter  = [0.0, 0.0, 0.0]; % isocenter of coil
    data.isTx           = 1; % if coil is also transmit (transceiver coil)
    data.parallelAP     = 0; % max parallelism in AP direction (0 is no parallel)
    data.parallelRL     = 0; % max parallelism in RL direction (0 is no parallel)
    data.parallelFH     = 0; % max parallelism in FH direction (0 is no parallel)
    % coil does not have data, treat it as optimal
    data.numCoils       = 1;
    data.coilConfig     = 'optimal';
    data.b1mIsocenter   = 1.0;
    data.b1pIsocenter   = 1.0;
    defaultCoil         = data.name;
    defaultIndex        = 1;
    b1mScaling          = data.b1mIsocenter;
    b1pScaling          = data.b1pIsocenter;

    % maps data
    maps.dim        = [];  % dimensions nx, ny, nz
    maps.spatial    = [];  % spatial coordinates
    maps.b1mSOS     = [];  % b1m Sum-of-Squares map
    maps.b1pSOS     = [];  % b1p Sum-of-Squares map
    maps.bx         = [];  % bx component maps
    maps.by         = [];  % by component maps
    
    % sar data and map
    sar.estSAR = 'N/A';
    sar.avgSAR = 'N/A';
    sar.cm3SAR = 'N/A';
    sar.g10SAR = 'N/A';
    sar.voxSAR = 'N/A';
    sar.mapSAR = [];

    % assign to struct
    coilSystem.coilModel{cc}.data = data;
    coilSystem.coilModel{cc}.maps = maps;
    coilSystem.coilModel{cc}.sar  = sar;
    
    fprintf(fid, '\n  %20s coil loaded -- parallelism AP=%d/RL=%d/FH=%d',...
        data.name, data.parallelAP, data.parallelRL, data.parallelFH);
    
end

%% prepare coil system struct with basic info common to all coils
coilSystem.numModels    = numModels;
coilSystem.isocenter    = systemIsocenter; % system isocenter
coilSystem.b1mScaling   = b1mScaling; % ref. scaling of the b1m maps
coilSystem.b1pScaling   = b1pScaling; % ref. scaling of the b1p maps
% active coils
coilSystem.activeTx     = defaultCoil; % active transmit name
coilSystem.indexTx      = defaultIndex; % active transmit index
coilSystem.activeRx     = defaultCoil; % active receive name
coilSystem.indexRx      = defaultIndex; % active receive index
% default transmit coil
coilSystem.defaultTx    = defaultCoil; % default transmit
coilSystem.defaultIndex = defaultIndex; % default transmit index
% transmit and receive modes
coilSystem.txMode       = 'ideal'; % mode: ideal / b1pSOS / pTX
coilSystem.rxMode       = 'pRX'; % mode: ideal / b1mSOS / pRX

%% report
tTotal = toc(tTotal);
fprintf(fid, '\n%s : %d coil models loaded into coilSystem',...
    functionName, coilSystem.numModels );
fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
fprintf(fid, '\n');
