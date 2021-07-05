function [acqData] = initializeAcqData()
%
% DATA.ACQUISITIONDATA.INITIALIZEACQDATA
%
%	Function that initializes basic acquisition data structure.
%   Returns a acqData with default fields.
%
%
% INPUT
%   None
%
% OUTPUT
%   acqData   structure acquisition data
%
%========================  CORSMED AB © 2020 ==============================
%

% basic info of sequence
acqData.pulseSeqSchem       = 0;
acqData.pulseSeqNum         = 0;
acqData.pulseSeqFamilyName  = '';
acqData.pulseSeqType        = 0;
acqData.pulseSeqAnalytical  = 0;
acqData.isCartesian         = 1;
acqData.is3D                = 0;
    
% basic info of acquisition
acqData.gamma           = 42.577478518e6;
acqData.matrixX         = 128;
acqData.matrixY         = 128;
acqData.matrixZ         = 1;
acqData.fovFE           = 0.200;
acqData.fovPE           = 0.200;
acqData.fovSE           = 0.200;
acqData.numFE           = 128;
acqData.numPE           = 128;
acqData.numSE           = 128;

acqData.sliceThickness  = 3e-3;
acqData.sliceGap        = 0.0;
acqData.numSlices       = 1;
acqData.zSliceCoord     = 0;

acqData.TE              = 10e-7; % echo time
acqData.TR              = 4.000; % local TR, rep duratioon
acqData.shotTR          = 4.000; % TR of shot, for example for MR-RAGE

acqData.rxBW         	= 250e3;
acqData.pixelBW       	= 781.25;

acqData.foldoverSuppr 	= 'no';
acqData.foldoverDir    	= 'AP';

% sampling factors
%  FE factor is used to mimic a LPF, and it is fixed
%  PE factor is either 1 or 2 depending on the foldoverSuppresion
%  SE factor is not used (fixed to 1)
acqData.samplingFactorFE    = 2; % BW multiplier factor for readout
acqData.samplingFactorPE    = 1; % foldover suppresion factor for phase enc
acqData.samplingFactorSE    = 1; % foldover suppresion factor for slice enc

acqData.NEX             = 1;

acqData.reconstructor   = 'coremri_fft';
acqData.parallelImaging = 'no';
acqData.rFactor         = 2;
acqData.partialFourier  = 'no';
acqData.fFactor         = 0.6;

acqData.numEchoes       = 1;
acqData.deltaTimeOOP    = 2.3e-3;
acqData.kspaceshift     = 'no';
acqData.ETL             = 1;
acqData.concatenations  = 0;
acqData.dummySSFP       = 10;
acqData.prepBSSFP       = 'ramp';
acqData.fatsat          = 'none';
acqData.interleaveSlice = 0; % apply or not interleaving for the slices

acqData.spoiled         = 0;
acqData.forceMinTE      = 0;
acqData.forceMinTR      = 0;
acqData.usePhaseSS      = 0; % use Phase for slice selection, otherwise freq
acqData.encOrder        = 'default'; % encoding order: default / reversed / centric / reversedCentric / ascending / descending

acqData.xPadFactor      = 1; % K-space padding factor before IFT, to increase image resolution
acqData.yPadFactor      = 1;
acqData.zPadFactor      = 1;

% plane info
acqData.pointsAll           = [];
acqData.pointsAllFrontend   = [];

acqData.vps             = 1;