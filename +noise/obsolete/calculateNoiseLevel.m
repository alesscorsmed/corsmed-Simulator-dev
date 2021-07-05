function [acquisition] = calculateNoiseLevel( ...
    acquisition, expControl, mrSystem )
%
% NOISE.CALCULATENOISELEVEL
%
%	generates noise data.
%
% INPUT
%   experiment          struct with experiment info from frontend
%   sessionData         solution struct with initial data
%
% OUTPUT
%   expControl          initialized experiment control struct
%   acquisition         initialized acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'noise.calculateNoiseLevel';
if (nargin < 3)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% get noise Level from DB
sqlQuery = ['SELECT selected_value FROM ',...
    'edt_tool_local.global_configuration WHERE',...
    ' name=''noise'''];
sqlQueryResults = exec(expControl.connLocalDB, sqlQuery);
noiseInfo   = fetch(sqlQueryResults);
applyNoise  = str2num(noiseInfo.Data{1,1});

%% get sequence related reference noise parameters
sqlQuery = ['SELECT voxelsizexref, voxelsizeyref, voxelsizezref, bwref,',...
    ' b0ref, nexref, kspacexref, kspaceyref, kspaceyref, noiseLevel',...
    ' FROM edt_tool_local.noise_level WHERE pulse_seq_scheme =',...
    num2str(acquisition.data.pulseSeqSchem)];
sqlQueryResults = exec(expControl.connLocalDB, sqlQuery);
refNoiseInfo    = fetch(sqlQueryResults);
refVoxelSizeX   = str2double(refNoiseInfo.Data{1,1});
refVoxelSizeY   = str2double(refNoiseInfo.Data{1,2});
refVoxelSizeZ   = str2double(refNoiseInfo.Data{1,3});
refBW           = str2double(refNoiseInfo.Data{1,4});
refB0           = str2double(refNoiseInfo.Data{1,5});
refNEX          = refNoiseInfo.Data{1,6};
refSizeX        = refNoiseInfo.Data{1,7};
refSizeY        = refNoiseInfo.Data{1,8};
refSizeZ        = refNoiseInfo.Data{1,9}; % Number of encoding slices
noiseLevel      = str2double(refNoiseInfo.Data{1,10}); % actual noise level
refVoxelVolume  = refVoxelSizeX*refVoxelSizeY*refVoxelSizeZ;
refEncSize      = refSizeX*refSizeY*refSizeZ;

if ~applyNoise
    noiseLevel = 0.0;
end

%% get relevant acquisition parameters
sizeX       = acquisition.data.numFE;
sizeY       = acquisition.data.numPE;
voxelSizeX  = acquisition.data.fovFE / sizeX;
voxelSizeY  = acquisition.data.fovPE / sizeY;
if acquisition.data.is3D
    sizeZ       = acquisition.data.numSlices;
    voxelSizeZ  = ( sizeZ * acquisition.data.sliceThickness ... 
        + (sizeZ-1) * acquisition.data.sliceGap )/ sizeZ;
else
    sizeZ       = 1;
    voxelSizeZ  = acquisition.data.sliceThickness;
end
fovVolume   = acquisition.data.fovFE*acquisition.data.fovPE*acquisition.data.fovSE;
voxelVolume = voxelSizeX*voxelSizeY*voxelSizeZ;
encSize     = sizeX*sizeY*sizeZ;
NEX         = acquisition.data.NEX;
BW          = acquisition.data.rxBW;
B0          = mrSystem.b0;

%% relative factors
factorBW    = BW/refBW;
factorB0    = B0/refB0;
factorVOL   = voxelVolume/refVoxelVolume;
factorNEX   = NEX/refNEX;
factorSize  = sizeX*sizeY*sizeZ/(refSizeX*refSizeY*refSizeZ);

%% compute the std noise
% % original function was
%
% noiseStd = anatomicalModelSlices*noiseLevel*sqrt(str2double(struct_exper.BW)/bwref)/...
%     ((B0/b0ref)*(voxelVolume/voxelVolumeRef)*...
%     sqrt((NEX/nexref)*((kspace(1,1)*kspace(1,2)*kspace(1,3))/(kspacexref*kspaceyref*kspacezref))));
%
% but I do not understand the anatomicalModelSlices factor,
% with anatomicalModelSlices the number of z isochromats
% of the model to simulate
%
noiseSD = 1000*noiseLevel/(factorB0*factorVOL) * ...
    sqrt( factorBW / (factorNEX*factorSize) );
% NOTE: 1000 scale is because in the old version the sliceThickness
% is multiplied by 1e-3 (line 32):
%   sliceThickness  = str2double(struct_exper.sliceThickness)*0.001; % from mm to m
% Which is weird because sliceThickness comes in m

%% reference SNR
refSNR = (refB0*refVoxelVolume/noiseLevel) * ...
    sqrt( refEncSize*refNEX/refBW );

%% current SNR (corrected by acceleration factors)
accEncSize = encSize;
if ~strcmpi(acquisition.data.parallelImaging,'no')
    accEncSize = accEncSize/acquisition.data.rFactor;
else
    if ~strcmpi(acquisition.data.partialFourier,'no')
        accEncSize = accEncSize*acquisition.data.fFactor;
    end
end
SNR = (B0*voxelVolume/noiseLevel) * ...
    sqrt( accEncSize*NEX/BW );

%% store data and return
acquisition.noise.noiseSD       = noiseSD;
acquisition.noise.relativeSNR   = round((SNR/refSNR)*100);
acquisition.noise.SNR           = SNR;
acquisition.noise.noiseLevel    = noiseLevel;
acquisition.noise.fovVolume     = fovVolume;
acquisition.noise.voxelVolume   = voxelVolume;
acquisition.noise.encSize       = encSize;
acquisition.noise.accEncSize    = accEncSize;
acquisition.noise.NEX           = NEX;
acquisition.noise.BW            = BW;
acquisition.noise.B0            = B0;
acquisition.noise.refSNR        = refSNR;
acquisition.noise.refVoxelVolume= voxelVolume;
acquisition.noise.refEncSize    = refEncSize;
acquisition.noise.refNEX        = refNEX;
acquisition.noise.refBW         = refBW;
acquisition.noise.refB0         = refB0;

