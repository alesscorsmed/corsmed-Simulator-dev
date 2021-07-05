function [pulseSequence] = se( acqData, rfData, mrSystem, dt, peEncodings )
%
% SPINTWIN.TEST.SEQ.SE
%
%	Generates a Spin Echo encoding:
%       RF90 + GR + RF180 + GR
%
% INPUT
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:test:seq:se';
% args
if (nargin<1) || isempty(acqData)
    acqData.rxBW                = 200e3;  % in Hz (this is a feature of the hardware)
    acqData.fovFE               = 0.3;
    acqData.fovPE               = 0.3;
    acqData.fovSE               = 0.3;
    acqData.numFE               = 128;
    acqData.numPE               = 128;
    acqData.numSE               = 128;
    acqData.samplingFactorFE    = 2;
    acqData.samplingFactorPE    = 1;
    acqData.samplingFactorSE    = 1;    
    acqData.echoTime            = 15e-3;
    acqData.repTime             = 4.000;
    acqData.preEncoding         = 1; % applies FE RW before 180
end
if (nargin<2) || isempty(rfData)
    rfData.type             = 'sinc';
    rfData.flipAngle      	= 90; % in degrees
    rfData.phase         	= 0.0; % in rads
    rfData.duration         = 1e-3;
    rfData.cycles           = 2;
    rfData.sliceThickness   = 3e-3;
    rfData.doSliceSelect 	= 0;
    rfData.preRWScale    	= 0;
    rfData.postRWScale    	= 0;
    rfData.keepTimingSS     = 0;
end
if (nargin<3) || isempty(mrSystem)
    mrSystem.b0             = 1.5000;
    mrSystem.maxGStrenght   = 0.0400;
    mrSystem.SlewRate       = 0.0;
end
if (nargin<4) || isempty(dt)
    dt = 1e-6;
end
if (nargin<5) || isempty(peEncodings)
    peEncodings = 0;
end

% constant
gamma = 42.577478518e6;

if acqData.preEncoding
    
    %% generate the encoding with readout
    [encTime,signalFE,signalPE,signalSE,signalRX,areaFE,dtEnc,effTE,rxLimits] = ...
        sequence.encodings.cartesianPhaseBalanced( ...
        acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
        acqData.numFE, acqData.numPE, acqData.numSE,...
        acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
        mrSystem.maxGStrenght, mrSystem.SlewRate, gamma, dt );
    
    %% generate Pre-Encoding Gradient
    [rwTime,rwSignal,~,~,~] = sequence.waveforms.grTrapArea(...
        areaFE/2, mrSystem.maxGStrenght, mrSystem.SlewRate, dt );
    
else
    
    %% generate a semi-balanced encoding (with FE RW)
    [encTime,signalFE,signalPE,signalSE,signalRX,dtEnc,effTE,rxLimits] = ...
        sequence.encodings.cartesianSemiBalanced( ...
        acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
        acqData.numFE, acqData.numPE, acqData.numSE,...
        acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
        mrSystem.maxGStrenght, mrSystem.SlewRate, gamma, dt );
    
    rwTime   = 10e-6;
    rwSignal = 0.0;
    
end

%% generate the waveforms for the RF 90 (use main RF info)
rfData.flipAngle = 90;
rfData.phase     = -pi/2;
[rf90Time,rfm90,rfp90,rff90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( rfData,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    gamma, dt );

%% generate the waveforms for the RF 180 (use refocusing RF info)
rfData.flipAngle = 180;
rfData.phase     = 0.0;
[rf180Time,rfm180,rfp180,rff180,grSlice180,rfLimits180] = ...
    sequence.excitations.sliceSelectRF( rfData,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    gamma, dt );

%% verify feasibility of TE and correctness of timing
minTE90  = rf90Time(end)/2 + rwTime(end) + rf180Time(end)/2;
tWait90  = acqData.echoTime - minTE90;
minTE180 = rf180Time(end)/2+ effTE;
tWait180 = acqData.echoTime - minTE180;
minTE = max(minTE180,minTE90);
if tWait90 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        acqData.echoTime*1e3, minTE*1e3 );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
if tWait180 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        acqData.echoTime*1e3, minTE*1e3 );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
      
%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRF180   = length(rf180Time);
numTimesRW      = length(rwTime);
numTimesRO      = length(encTime);
numTimesWait90  = round(tWait90/dt);
numTimesWait180 = round(tWait180/dt);

% start/end of RF90
iniRF90         = 1;
endRF90         = iniRF90 + numTimesRF90-1;
rfLimits90      = rfLimits90 + iniRF90-1;
% start/end of 1st gradient
iniRwGr         = endRF90 + 1;
endRwGr         = iniRwGr+numTimesRW-1;
% start/end of RF180
iniRF180        = endRwGr + numTimesWait90 + 1;
endRF180        = iniRF180 + numTimesRF180-1;
rfLimits180     = rfLimits180 + iniRF180-1;
% start/end of readout gradient
iniRoGr         = endRF180 + numTimesWait180 + 1;
endRoGr         = iniRoGr+numTimesRO-1;
% total time of sequence
numSteps        = endRoGr + 1;

%% pulse sequence
pulseSequence = data.simulation.initializeSequence();
pulseSequence.type         = 'SE';
pulseSequence.gamma        = gamma;

% waveforms
repetition.time      = zeros(numSteps,1);
repetition.rxSignal  = zeros(numSteps,1); % receiver readout
repetition.swcSignal = zeros(numSteps,1); % software crusher
repetition.gxSignal  = zeros(numSteps,1); % x gradient
repetition.gySignal  = zeros(numSteps,1); % y gradient
repetition.gzSignal  = zeros(numSteps,1); % z gradient
repetition.rfmSignal = zeros(numSteps,1); % RF magnitude
repetition.rfpSignal = zeros(numSteps,1); % RF phase
repetition.rffSignal = zeros(numSteps,1); % RF frequency
repetition.gdwSignal = zeros(numSteps,3); % diffusion gradients

%% assign blocks
repetition.time(:)     = dt*(1:numSteps);
%% RF
repetition.rfmSignal(iniRF90:endRF90)    = rfm90;
repetition.rfpSignal(iniRF90:endRF90)    = rfp90;
repetition.rffSignal(iniRF90:endRF90)    = rff90;
repetition.gzSignal(iniRF90:endRF90)     = grSlice90;
%% Refocusing
repetition.rfmSignal(iniRF180:endRF180)  = rfm180;
repetition.rfpSignal(iniRF180:endRF180)  = rfp180;
repetition.rffSignal(iniRF180:endRF180)  = rff180;
repetition.gzSignal(iniRF180:endRF180)   = grSlice180;
%% gradients RW
repetition.gxSignal(iniRwGr:endRwGr,1) =...
    repetition.gxSignal(iniRwGr:endRwGr,1) + rwSignal(:);
%% gradients RO
repetition.gxSignal(iniRoGr:endRoGr,1) =...
    repetition.gxSignal(iniRoGr:endRoGr,1) + signalFE(:);
repetition.gySignal(iniRoGr:endRoGr,1) =...
    repetition.gySignal(iniRoGr:endRoGr,1) + signalPE(:);
repetition.gzSignal(iniRoGr:endRoGr,1) =...
    repetition.gzSignal(iniRoGr:endRoGr,1) + 0*signalSE(:);
%% Samples
numSamples = numel(signalRX);
repetition.rxSignal(iniRoGr:endRoGr) = 1:numSamples; % RO during whole encoding
rxLimits = [iniRoGr, endRoGr];
%% SWC
repetition.swcSignal(end) = 1;

%% for splitting into parts (RF vs Non-RF)
numParts = 1;
% 90 
if rfLimits90(1) > 1
    numParts = numParts + 1;
    partLimits = [1, rfLimits90(1)-1; rfLimits90];
    partType{1} = 'GR';
    partType{2} = 'RF';
else
    partLimits = rfLimits90;
    partType{1} = 'RF';
end
% first gradient
numParts = numParts + 1;
partLimits = [partLimits; rfLimits90(2)+1, rfLimits180(1)-1];
partType{numParts} = 'GR';
% 180
numParts = numParts + 1;
partLimits = [partLimits; rfLimits180];
partType{numParts} = 'RF';
% encoding with readouts
numParts = numParts + 1;
partLimits = [partLimits; rfLimits180(2)+1, numSteps];
partType{numParts} = 'RO';

% apply reps
numReps = numel(peEncodings);

%% first rep
pulseSequence.time      = repetition.time;
pulseSequence.rxSignal  = repetition.rxSignal;
pulseSequence.swcSignal = repetition.swcSignal;
pulseSequence.gxSignal  = repetition.gxSignal;
pulseSequence.gySignal  = repetition.gySignal * peEncodings(1);
pulseSequence.gzSignal  = repetition.gzSignal;
pulseSequence.rfmSignal = repetition.rfmSignal;
pulseSequence.rfpSignal = repetition.rfpSignal;
pulseSequence.rffSignal = repetition.rffSignal;
pulseSequence.gdwSignal = repetition.gdwSignal;

% prepare for more parts
repRxLimits   = rxLimits;
repPartLimits = partLimits;
repPartType   = partType;
repNumParts   = numParts;

%% more repetitions
TR = acqData.repTime;
for ii = 2:numReps
       
    pulseSequence.time(end) = (ii-1)*TR;      
    pulseSequence.time      = [ pulseSequence.time; repetition.time + (ii-1)*TR ];

    pulseSequence.swcSignal = [ pulseSequence.swcSignal; repetition.swcSignal];
    pulseSequence.gxSignal  = [ pulseSequence.gxSignal; repetition.gxSignal];
    pulseSequence.gySignal  = [ pulseSequence.gySignal; repetition.gySignal * peEncodings(ii)];
    pulseSequence.gzSignal  = [ pulseSequence.gzSignal; repetition.gzSignal];
    pulseSequence.rfmSignal = [ pulseSequence.rfmSignal; repetition.rfmSignal];
    pulseSequence.rfpSignal = [ pulseSequence.rfpSignal; repetition.rfpSignal];
    pulseSequence.rffSignal = [ pulseSequence.rffSignal; repetition.rffSignal];
    pulseSequence.gdwSignal = [ pulseSequence.gdwSignal; repetition.gdwSignal];
    
    repRxSignal = repetition.rxSignal;
    repRxSignal(repRxSignal > 0) = repRxSignal(repRxSignal > 0) + (ii-1)*numSamples;
    pulseSequence.rxSignal  = [ pulseSequence.rxSignal; repRxSignal];
    rxLimits = [rxLimits; repRxLimits + (ii-1)*numSteps];
    
    partLimits = [partLimits; repPartLimits+(ii-1)*numSteps];
    partType(numParts+1:numParts + repNumParts) = repPartType;
    numParts = numParts + repNumParts;
end


%% final info
pulseSequence.timeDiff     = reshape([dt; diff(pulseSequence.time(:))],[],1);
pulseSequence.totalTime    = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps     = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs       = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.effTE        = acqData.echoTime;
pulseSequence.rxLimits     = rxLimits; % starting / end indexes of Readouts
% assign parts
pulseSequence.numParts   = numParts; % number of parts
pulseSequence.partType   = partType; % type of part
pulseSequence.partLimits = partLimits; % index start/end of part
