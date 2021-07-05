function [acqNoise] = initializeAcqNoise()
%
% DATA.ACQUISITIONDATA.INITIALIZEACQNOISE
%
%	Function that initializes basic acquisition noise structure.
%   Returns a acqNoise with default fields.
%
%
% INPUT
%   None
%
% OUTPUT
%   acqNoise   structure acquisition noise
%
%========================  CORSMED AB © 2020 ==============================
%

acqNoise.noiseLevel    = 0;
acqNoise.noiseSD       = 0;
acqNoise.relativeSNR   = 0;
acqNoise.SNR           = 0;
acqNoise.fovVolume     = 0;
acqNoise.voxelVolume   = 0;
acqNoise.encSize       = 0;
acqNoise.accEncSize    = 0;
acqNoise.NEX           = 0;
acqNoise.BW            = 0;
acqNoise.B0            = 0;
acqNoise.refSNR        = 0;
acqNoise.refVoxelVolume= 0;
acqNoise.refEncSize    = 0;
acqNoise.refNEX        = 0;
acqNoise.refBW         = 0;
acqNoise.refB0         = 0;
