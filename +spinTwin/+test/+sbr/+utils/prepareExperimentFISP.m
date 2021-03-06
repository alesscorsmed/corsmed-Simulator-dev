function [pulseSequence,simModel,motionModel,mrSystem] = ...
    prepareExperimentFISP( fovX, fovY, fovZ, dx, dy, dz)
%
% SPINTWIN.SBR.EXAMPLE.UTILS.PREPAREEXPERIMENT
%
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

%% setup simulation object and motion content
[simModel, motionModel] = spinTwin.test.sbr.utils.prepObject(fovX, fovY, fovZ, dx, dy, dz);


%% define system main limits
mrSystem.b0             = 1.5000;
mrSystem.maxGStrenght   = 0.0600;
mrSystem.SlewRate       = 0.0;


%% define sequence (IR-FISP train)
[expControl]   = data.experiment.initializeExpControl();
[acquisition]  = data.experiment.initializeAcquisition();

% set up parameters for the sequence
acquisition.data.pulseSeqActualName = 'ftrain';
acquisition.data.pulseSeqFamilyName = 'fisp-train';
acquisition.data.gamma              = 42.577478518e6;
acquisition.data.matrixX            = 64;
acquisition.data.matrixY            = 64;
acquisition.data.matrixZ            = 1;
acquisition.data.fovFE              = fovX;
acquisition.data.fovPE              = fovY;
acquisition.data.fovSE              = fovZ;
acquisition.data.numFE              = 64;
acquisition.data.numPE              = 64;
acquisition.data.numSE              = 1;
acquisition.data.zSliceCoord        = 0.0;
acquisition.data.sliceThickness     = fovZ;
acquisition.data.TE                 = 5e-3;
acquisition.data.TR                 = 10e-3;
acquisition.prepIR.TI               = 50e-3;

shots                               = acquisition.data.numPE * 4;
acquisition.data.rfmod              = abs(sin(2*(1:shots)'*pi/180));

%br = @(x)abs(sin(2*(x)*pi/180));
%plot(1:64*10,br(1:64*10))


% apply encoding
[encoding] = encoder.generateEncodingPlan( ...
    acquisition, expControl );

% generate the sequence
[pulseSequence] = sequence.generatePulseSequence( ...
         acquisition, encoding, mrSystem, expControl );


