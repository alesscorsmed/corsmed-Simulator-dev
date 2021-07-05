function [pulseSequence,simModel,motionModel,mrSystem] = ...
    prepareExperimentDyn( fovX, fovY, fovZ, dx, dy, dz)
%
% SPINTWIN.SBR.EXAMPLE.UTILS.PREPAREEXPERIMENT
%
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

%% prepare the data required

%% model and data
b0          = 1.0;
mu          = 1.0;
pd          = 1.0;
t2Limits    = [0.100,  1.000];
t1Limits    = [0.150,  3.000];
numT2       = 4;
numT1       = 4;

%% generate a compartment model
t1Value = mean(t1Limits);
t2Value = mean(t2Limits);

% generate spin Model
[spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
    fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );

% get and modify the simulation model
simModel = spinModel.slice{1}.model; clear spinModel;

% get dimensions
[nX, nY, nZ, ~] = size(simModel.r3D);
% prepare tissue properties
simModel.pd	= pd*ones(nX, nY, nZ);
simModel.r1	= zeros(nX, nY, nZ);
simModel.r2	= zeros(nX, nY, nZ);
simModel.bi	= zeros(nX, nY, nZ);
simModel.cs	= zeros(nX, nY, nZ);

% t1/t2 values
t2Values = logspace(log10(min(t2Limits)), log10(max(t2Limits)), numT2);
t1Values = logspace(log10(min(t1Limits)), log10(max(t1Limits)), numT1);

% assign t1/t2 pairs to each compartment
localNx = round(nX/numT1);
localNy = round(nY/numT2);
for ii=1:numT1
    idxII = 1+(ii-1)*localNx:min(ii*localNx,nX);
    for jj=1:numT2
        idxJJ = 1+(jj-1)*localNy:min(jj*localNy,nY);
        t2Val = t2Values(jj);
        t1Val = t1Values(ii);
        t1Val = max(t1Val, t2Val + 10e-3);
        simModel.r1(idxII,idxJJ,:) = 1/t1Val;
        simModel.r2(idxII,idxJJ,:) = 1/t2Val;
    end
end

%% mrSystem
mrSystem.b0             = 1.5000;
mrSystem.maxGStrenght   = 0.0600;
mrSystem.SlewRate       = 0.0;

%% motion
motionModel             = [];



%% sequence: IR SS FSE
[expControl]   = data.experiment.initializeExpControl();
[acquisition]  = data.experiment.initializeAcquisition();

% set up parameters for the sequence
acquisition.data.pulseSeqActualName = 'Joda';
acquisition.data.pulseSeqFamilyName = 'bssfp-train';
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
acquisition.data.TE                 = 7e-3;
acquisition.data.TR                 = 15e-3;
acquisition.prepIR.TI               = 50e-3;

% apply encoding
[encoding] = encoder.generateEncodingPlan( ...
    acquisition, expControl );

% generate the sequence
[pulseSequence] = sequence.generatePulseSequence( ...
         acquisition, encoding, mrSystem, expControl );


