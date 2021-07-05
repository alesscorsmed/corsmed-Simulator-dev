function [acquisition] = calculateNoiseLevel( acquisition, expControl, B0 )
%
% NOISE.CALCULATENOISELEVEL
%
%	generates noise data.
%
% INPUT
%   acquisition 	initialized acquisition structure
%   B0              main field strength in T
%
% OUTPUT
%   acquisition  	updated acquisition structure with noise
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'noise.calculateNoiseLevel';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end
if (nargin < 2) || isempty(B0)
    B0 = 1.5;
end

%% extract noise parameters
refBW           = acquisition.noise.refBW;
refB0           = acquisition.noise.refB0;
refNEX          = acquisition.noise.refNEX;
refVoxelVolume  = acquisition.noise.refVoxelVolume;
refEncSize      = acquisition.noise.refEncSize;

if expControl.simulation.applyNoise
    noiseLevel 	= acquisition.noise.noiseLevel;
else
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

%% relative factors
factorBW    = BW/refBW;
factorB0    = B0/refB0;
factorVOL   = voxelVolume/refVoxelVolume;
factorNEX   = NEX/refNEX;
factorSize  = encSize/(refEncSize);

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
acquisition.noise.fovVolume     = fovVolume;
acquisition.noise.voxelVolume   = voxelVolume;
acquisition.noise.encSize       = encSize;
acquisition.noise.accEncSize    = accEncSize;
acquisition.noise.NEX           = NEX;
acquisition.noise.BW            = BW;
acquisition.noise.B0            = B0;
acquisition.noise.refSNR        = refSNR;
acquisition.noise.refVoxelVolume= refVoxelVolume;
acquisition.noise.refEncSize    = refEncSize;
acquisition.noise.refNEX        = refNEX;
acquisition.noise.refBW         = refBW;
acquisition.noise.refB0         = refB0;
acquisition.noise.noiseLevel    = noiseLevel;

