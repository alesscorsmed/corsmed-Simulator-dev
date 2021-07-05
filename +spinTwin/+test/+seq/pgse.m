function [pulseSequence, grData] = pgse( grData, rfData, mrSystem, dt )
%
% SPINTWIN.TEST.SEQ.PGSE
%
%	Generates a Pulse Gradient Spin Echo encoding:
%       RF90 + GR + RF180 + GR
%
% INPUT
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:test:seq:pgse';
% args
if (nargin<1) || isempty(grData)
    grData.echoTime         = 50e-3;
    grData.duration      	= 25e-3;
    grData.amplitude     	= 0.040;
    grData.direction        = 'x'; % encoding direction 'x' / 'y' / 'z'
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

% constant
gamma = 42.577478518e6;

%% generate Pulsed Gradient
[pgTime,pgSignal,pgArea,~,pgPlatLimits] = sequence.waveforms.grTrapPlateau(...
        grData.duration, grData.amplitude, mrSystem.SlewRate, dt);
% compute remaining data    
grData.area     = pgArea; % gradient area
grData.time     = pgTime(end); % total time
grData.plateau  = diff(pgTime(pgPlatLimits)); % plateau time
grData.beta     = (2*pi*gamma)^2 *(grData.area^2) *(grData.echoTime-grData.time/3);

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
minTE90  = rf90Time(end)/2 + pgTime(end) + rf180Time(end)/2;
tWait90  = grData.echoTime - minTE90;
minTE180 = rf180Time(end)/2+ pgTime(end);
tWait180 = grData.echoTime - minTE180;
minTE = max(minTE180,minTE90);
if tWait90 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        grData.echoTime*1e3, minTE*1e3 );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
if tWait180 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        grData.echoTime*1e3, minTE*1e3 );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
      
%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRF180   = length(rf180Time);
numTimesGrad    = length(pgTime);
numTimesWait90  = round(tWait90/dt);
numTimesWait180 = round(tWait180/dt);
numTimesTAU     = round(grData.echoTime/dt);

% start/end of RF90
iniRF90         = 1;
endRF90         = iniRF90 + numTimesRF90-1;
rfLimits90      = rfLimits90 + iniRF90-1;
% start/end of 1st gradient
iniFwGr         = endRF90 + 1;
endFwGr         = iniFwGr+numTimesGrad-1;
% start/end of RF180
iniRF180        = endFwGr + numTimesWait90 + 1;
endRF180        = iniRF180 + numTimesRF180-1;
rfLimits180     = rfLimits180 + iniRF180-1;
% start/end of rewind gradient
iniRwGr         = endRF180 + 1;
endRwGr         = iniRwGr+numTimesGrad-1;
% total time of sequence
numSteps        = max( round(iniRF90 + numTimesRF90/2 + 2*numTimesTAU), ...
    endRF180 + numTimesWait180 );

%% pulse sequence
pulseSequence = data.simulation.initializeSequence();
pulseSequence.type         = 'PGSE encoding';
pulseSequence.gamma        = gamma;

% waveforms
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.timeDiff  = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1); % receiver readout
pulseSequence.swcSignal = zeros(numSteps,1); % software crusher
pulseSequence.gxSignal  = zeros(numSteps,1); % x gradient
pulseSequence.gySignal  = zeros(numSteps,1); % y gradient
pulseSequence.gzSignal  = zeros(numSteps,1); % z gradient
pulseSequence.rfmSignal = zeros(numSteps,1); % RF magnitude
pulseSequence.rfpSignal = zeros(numSteps,1); % RF phase
pulseSequence.rffSignal = zeros(numSteps,1); % RF frequency
pulseSequence.gdwSignal = zeros(numSteps,3); % diffusion gradients

%% assign blocks
pulseSequence.time(:)     = dt*(1:numSteps);
pulseSequence.timeDiff(:) = dt;
%% RF
pulseSequence.rfmSignal(iniRF90:endRF90)    = rfm90;
pulseSequence.rfpSignal(iniRF90:endRF90)    = rfp90;
pulseSequence.rffSignal(iniRF90:endRF90)    = rff90;
pulseSequence.gzSignal(iniRF90:endRF90)     = grSlice90;
%% Refocusing
pulseSequence.rfmSignal(iniRF180:endRF180)  = rfm180;
pulseSequence.rfpSignal(iniRF180:endRF180)  = rfp180;
pulseSequence.rffSignal(iniRF180:endRF180)  = rff180;
pulseSequence.gzSignal(iniRF180:endRF180)   = grSlice180;
%% PG gradients (Fw and Rw)
switch lower(grData.direction)
    case 'x'
        pulseSequence.gxSignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gxSignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gxSignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gxSignal(iniRwGr:endRwGr,1) + pgSignal(:);
    case 'y'
        pulseSequence.gySignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gySignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gySignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gySignal(iniRwGr:endRwGr,1) + pgSignal(:);
    case 'z'
        pulseSequence.gzSignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gzSignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gzSignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gzSignal(iniRwGr:endRwGr,1) + pgSignal(:);
    otherwise
        msg = sprintf('gradient direction %s not supported',...
            grData.direction );
        ME = MException(['error:',functionName], '%s', msg);
        throw(ME);
end
pulseSequence.rxSignal(end)    = 1; % single RO at end
rxLimits = [numSteps, numSteps];

%% pulse sequences info
pulseSequence.totalTime    = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps     = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs       = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.effTE        = grData.echoTime;
pulseSequence.rxLimits     = rxLimits; % starting / end indexes of Readouts

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
% second gradient
numParts = numParts + 1;
partLimits = [partLimits; rfLimits180(2)+1, endRwGr];
partType{numParts} = 'GR';
% rest may have readouts
if pulseSequence.numSteps > endRwGr
    numParts = numParts + 1;
    partLimits = [partLimits; endRwGr+1, pulseSequence.numSteps];
    partType{numParts} = 'RO';
end
% assign
pulseSequence.numParts   = numParts; % number of parts
pulseSequence.partType   = partType; % type of part
pulseSequence.partLimits = partLimits; % index start/end of part


