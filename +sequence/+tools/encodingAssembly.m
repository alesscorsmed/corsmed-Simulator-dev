function [time, signalFE, signalPE, signalSE, signalRX, rxLimits] = ...
    encodingAssembly( encTime, encSignal, rxSignal, rxLimits, ...
    feRiseSteps, feRwSignal, peRwSignal, seRwSignal, ...
    feRwScale, peRwScale, seRwScale, dt, expControl )
%
% SEQUENCE.TOOLS.ENCODINGASSEMBLY
%
%	Assemble multiple echoes into a single encoding
%
% INPUT
%   numRX           number of samples to place per interval
%   dtBW            time step for the sampling
%   time            time points
%   signalFE        frequency encoding gradient signal
%   platLimits      limits for the placement (plateau)
%   expControl      
%
% OUTPUT
%   newTime         time points, starting in dt
%   newSignalFE     frequency encoding gradient signal
%   newPlatLimits   limits for the placement (plateau)
%   newSignalRX     readout signal
%   newRxLimits     start and end indexes of the first readout
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:tools:encodingAssembly';
% check args
if (nargin < 12)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 13) || isempty(expControl)
    expControl.debug.debugMode = 0;
    expControl.debug.debugFile = [];
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

% get dimensions of each signals
lenPeRw     = numel(peRwSignal);
lenSeRw     = numel(seRwSignal);
lenPhase    = max(lenPeRw,lenSeRw);
lenFeRw     = numel(feRwSignal);
lenFeRise   = feRiseSteps;
lenFeSignal = numel(encTime);
% number of total steps
lenEnc = lenFeSignal + nnz(feRwScale)*lenFeRw; % total steps for FE
% pre encoding delay (number of extra steps, accounting for overlap)
if seRwScale(1) || peRwScale(1)
    lenFePreDL = lenPhase - (lenFeRise-1);
    if feRwScale(1)
        lenFePreDL = lenFePreDL - lenFeRw;
    end
    lenFePreDL = max(0,lenFePreDL);
else
    lenFePreDL = 0;
end
% post encoding delay (number of extra steps, accounting for overlap)
if seRwScale(2) || peRwScale(2)
    lenFePostDL = lenPhase - (lenFeRise-1);
    if feRwScale(2)
        lenFePostDL = lenFePostDL - lenFeRw;
    end
    lenFePostDL = max(0,lenFePostDL);
else
    lenFePostDL = 0;
end
% total number of points
numTimePoints = lenEnc + lenFePreDL + lenFePostDL;  

%% combine into signals
% create signals
signalFE = zeros(numTimePoints,1);
signalPE = zeros(numTimePoints,1);
signalSE = zeros(numTimePoints,1);
signalRX = zeros(numTimePoints,1);
time     = zeros(numTimePoints,1);

if seRwScale(1) || peRwScale(1)
    % populate phase rewinders at starting point
    signalPE(1:lenPeRw)   	= peRwScale(1)*peRwSignal(:);
    signalSE(1:lenSeRw)   	= seRwScale(1)*seRwSignal(:);
end

idxShift  =  lenFePreDL; % update the shift to current position
if lenFeRw > 0
    % pre FE rewinder
    idxTarget               = (1:lenFeRw) + idxShift;
    signalFE(idxTarget)     =  feRwScale(1)*feRwSignal(:);
    idxShift                =  idxShift + lenFeRw; % update shift
end

if idxShift > 0
    % common time step
    time(1:idxShift)        =  dt*(1:idxShift); % common time step is dt
    timeShift               =  time(idxShift);
else
    timeShift               = 0;
end

% FE encoding
if lenFeSignal > 0
    idxTarget               = (1:lenFeSignal) + idxShift;
    signalFE(idxTarget)     =  encSignal(:);
    signalRX(idxTarget)     =  rxSignal(:);
    time(idxTarget)         =  encTime + timeShift;
    rxLimits                =  rxLimits + idxShift; % shift the rxLimits indexes
    idxShift                =  idxShift + lenFeSignal; % update shift
    timeShift               =  time(idxShift);
end

% remaining time is in dt steps
time(idxShift+1:end)    = timeShift + dt*(1:numTimePoints-idxShift);

% post rewinders (if needed)
if feRwScale(2)
    % post FE rewinder
    idxTarget               = (1:lenFeRw) + idxShift;
    signalFE(idxTarget)     = feRwScale(2)*feRwSignal(:);
end
if seRwScale(2) || peRwScale(2)
    % post phase rewinders
    signalSE(numTimePoints+1 -(1:lenSeRw))	= seRwScale(2)*seRwSignal(:); % reverse from end
    signalPE(numTimePoints+1 -(1:lenPeRw))	= peRwScale(2)*peRwSignal(:); % reverse from end
end


%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end


