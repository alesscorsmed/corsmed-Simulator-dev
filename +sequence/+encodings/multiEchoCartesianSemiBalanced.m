function [time,signalFE,signalPE,signalSE,signalRX,dt,effTE,rxLimits] = ...
    multiEchoCartesianSemiBalanced( echoTimes, receiverBW, ...
    fovFE, fovPE, fovSE, numFE, numPE, numSE, ...
    samplingFactorFE, samplingFactorPE, samplingFactorSE, ...
    maxGStrenght, gradSlewrate, gamma, dt, expControl)
%
% SEQUENCE.ENCODINGS.MULTIECHOCARTESIANSEMIBALANCED
%
%	Generates a cartesian Semi-Balanced encoding (for example for GRE).
%   FE:  -HalfRewind  +Encoding -Encoding +Encoding ...
%   PE:   HalfRewind
%   SE:   HalfRewind
%
% INPUT
%   echoTimes       array with times between centers of Encodings
%   receiverBW      bandwidth of HW receiver in Hz (will determine dt)
%   fovFE           Field of View in the frequency encoding direction
%   fovPE           Field of View in the phase encoding direction
%   fovSE           Field of View in the 3D (slice) encoding direction
%   numFE           number of frequency encodings (samples in readout)
%   numPE           number of phase encodings
%   numSE           number of 3D encodings
%   maxGStrenght    max gradient amplitude in T/m
%   gradSlewrate    slew rate in T/m/s
%   gamma           gyromagnetic ratio
%   dt              max allowed dt
%   debugMode       flag for debugging
%   debugFile       name of file where to dump debugging info
%
% OUTPUT
%   time            time points, starting in dt
%   signalFE        frequency encoding gradient signal
%   signalPE        phase encoding gradient signal
%   signalSE        3D (slice) encoding gradient signal
%   signalRX        readout signal
%   dt              time step used
%   effTE           effective Time Echo: time to the center of the encoding
%   rxLimits        start and end indexes of the first readout
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.encodings.multiEchoCartesianSemiBalanced';

if (nargin < 1 || isempty(echoTimes))
    echoTimes = 2.3e-3;  % s
end
if (nargin < 2 || isempty(receiverBW))
    receiverBW = 200e3;  % in Hz (this is a feature of the hardware)
end
if (nargin < 3 || isempty(fovFE))
    fovFE   = 0.3;
end
if (nargin < 4 || isempty(fovPE))
    fovPE 	= 0.3;
end
if (nargin < 5 || isempty(fovSE))
    fovSE 	= 0.3;
end
if (nargin < 6 || isempty(numFE))
    numFE	= 256;
end
if (nargin < 7 || isempty(numPE))
    numPE	= 128;
end
if (nargin < 8 || isempty(numSE))
    numPE	= 128;
end
if (nargin < 9 || isempty(samplingFactorFE))
    samplingFactorFE = 1;
end
if (nargin < 10 || isempty(samplingFactorPE))
    samplingFactorPE = 1;
end
if (nargin < 11 || isempty(samplingFactorSE))
    samplingFactorSE = 1;
end
if (nargin < 12 || isempty(maxGStrenght))
    maxGStrenght = 0.030;
end
if (nargin < 13 || isempty(gradSlewrate))
    gradSlewrate = 150.0;
end
if (nargin < 14 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end
if (nargin < 15 || isempty(dt))
    dt = 1e-6; % max allowed dt
end
if (nargin < 16 || isempty(expControl))
    expControl.connLocalDB = [];
    expControl.application = 'unknown';
    expControl.debug.debugMode = 1;
    expControl.debug.debugFile = '';
end

% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% apply the sampling factor to each direction
% scale FOV and number of encodings in Freq direction: mimic of LPF
numFE   = samplingFactorFE*numFE;
fovFE   = samplingFactorFE*fovFE;
% scale bandwidth accordingly to increase sampling rate for simulation
simRxBW = samplingFactorFE*receiverBW;
% scale FOV and number of encodings in Phase direction: foldover related
numPE   = samplingFactorPE*numPE;
fovPE   = samplingFactorPE*fovPE;
% scale FOV and number of encodings in Slice direction: useless, always 1
numSE   = samplingFactorSE*numSE;
fovSE   = samplingFactorSE*fovSE;

%% time discretization as multiple of simulation receiverBW cycle and less than dt
dtBW = 1/simRxBW;
samplingRate = ceil(dtBW/dt);
dt = 1/(simRxBW*samplingRate);

%% generate FE gradient signal for Readout
gradAmp  = simRxBW/(gamma*fovFE);
% check if maximum gradient limit is trespassed in Readout
if gradAmp > maxGStrenght
    msg = sprintf( ['Gradient amplitude in readout (%.3fmT/m) '... 
        'exceeds system maximum (%.3fmT/m). '...
        'Consider reducing the BW.'],...
        gradAmp*1e3, maxGStrenght*1e3 );
    eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
end
% generate waveform
gradTime = numFE/simRxBW;
[feTime,feSignal,feArea,~,fePlatLimits] = ...
    sequence.waveforms.grTrapPlateau(gradTime,gradAmp,gradSlewrate,dt);

%% verify time between readouts
minEchoTime = feTime(end);
if minEchoTime > min(echoTimes)
    msg = sprintf( ['Times between Echoes (%.3fms) '...
        'are too short (minimum %.3fms). '...
        'Consider increasing time between echoes.'],...
        min(echoTimes)*1e3, minEchoTime*1e3 );
    eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
end

%% generate FE rewinder
[feRwTime,feRwSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(feArea/2,maxGStrenght,gradSlewrate,dt);

%% generate half of PE encoding
peArea = numPE/(gamma*fovPE)/2;
[peTime,peSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(peArea,maxGStrenght,gradSlewrate,dt);

%% generate half of 3D encoding
seArea = numSE/(gamma*fovSE)/2;
[seTime,seSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(seArea,maxGStrenght,gradSlewrate,dt);

%% construct single part of the Encoding
numEchoes   = length(echoTimes) + 1;
lenEchoDelay= round((echoTimes - minEchoTime)/dt);
lenPe       = length(peTime);
lenSe       = length(seTime);
lenPhase    = max(lenPe,lenSe);
lenFeRw     = length(feRwTime);
lenFeRise   = fePlatLimits(1);
lenFeSignal = length(feTime);
lenFeDelay  = lenPhase - (lenFeRise+lenFeRw); % steps need to add to FE
lenEnc      = lenFeRw + numEchoes*lenFeSignal + sum(lenEchoDelay(:)); % total tsteps
if lenFeDelay > 0
    % add extra time steps
    lenEnc = lenEnc + lenFeDelay; 
else
    lenFeDelay = 0;
end

%% combine parts into signals
% last active entry of the frequency encoding
numTimePoints   = lenEnc+1;
% create signals
signalFE = zeros(numTimePoints,1);
signalPE = zeros(numTimePoints,1);
signalSE = zeros(numTimePoints,1);
signalRX = zeros(numTimePoints,1);
time     = zeros(numTimePoints,1);
% populate
signalPE(1:lenPe)                   =  peSignal(:);
signalSE(1:lenSe)                   =  seSignal(:);
idxShift = lenFeDelay;
signalFE( (1:lenFeRw) + idxShift)   = -feRwSignal(:);
idxShift = idxShift + lenFeRw;
signalFE( (1:lenFeSignal) + idxShift)   =  feSignal(:);
% readout placement
echoFePlatLimits = fePlatLimits + idxShift;
roIdx = echoFePlatLimits(1):samplingRate:echoFePlatLimits(1)+samplingRate*numFE;
roIdx = roIdx(1:numFE);
roIdx = roIdx + floor( (echoFePlatLimits(2) - roIdx(end))/2 );
signalRX(roIdx) = (1:numFE);
% index limits of first readout
rxLimits = [roIdx(1), roIdx(end)];
% increase index shift
idxShift = idxShift + lenFeSignal;
% rest of echoes
lenEchoDelay(numEchoes) = 0;
for ii = 2:numEchoes
    % add delay between echoes
    idxShift = idxShift + lenEchoDelay(ii);
    % add FE for echo, with sign
    sign = (-1)^(ii-1);
    signalFE( (1:lenFeSignal) + idxShift) = sign*feSignal(:);
    % correct readout placement
    echoFePlatLimits = fePlatLimits + idxShift;
    % readout placement
    roIdx = echoFePlatLimits(1):samplingRate:echoFePlatLimits(1)+samplingRate*numFE;
    roIdx = roIdx(1:numFE);
    roIdx = roIdx + floor( (echoFePlatLimits(2) - roIdx(end))/2 );
    signalRX(roIdx) = (1:numFE) + (ii-1)*numFE;
    % index limits of first readout
    rxLimits = [rxLimits; roIdx(1), roIdx(end)];
    % increase shift
    idxShift = idxShift + lenFeSignal;
end
% time
time(1:end) = dt*(1:numTimePoints);
% effective Echo timing
effTE  = mean(time(rxLimits(1,:)));

% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fms', 1e3*time(end));
    fprintf(fid, '\n  Eff. Time Echo      %.3fms', effTE*1e3 );
    fprintf(fid, '\n  Rx Bandwidth        %.3fKHz', simRxBW*1e-3);
    fprintf(fid, '\n  Rx Sampling Rate    %.3fus', 1e6/simRxBW);
    fprintf(fid, '\n  Time step           %.3fus', dt*1e6);
    fprintf(fid, '\n  Number steps        %d', numTimePoints);
    fprintf(fid, '\n  Number readouts     %d', nnz(signalRX));
    fprintf(fid, '\n  Number Echoes       %d', numEchoes);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
%         % plot
%         plot(time, signalFE)
%         hold on
%         plot(time, signalPE)
%         hold on
%         plot(time, signalSE)
%         hold on
%         plot(time(signalRX>0), signalFE(signalRX>0), 'o')
    
end

