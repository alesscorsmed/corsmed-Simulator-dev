function [time, signalFE, signalRX, rxLimits ] = ...
    multiEchoAssembly( feTime, feSignal, feRxSignal, feRxLimits,...
    echoTimes, reversePolarity, feFullRwTime, feFullRwSignal, expControl )
%
% SEQUENCE.TOOLS.MULTIECHOASSEMBLY
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
functionName = 'sequence:tools:multiEchoAssembly';
% check args
if (nargin < 4) ...
        || isempty(feTime) || isempty(feSignal) ...
        || isempty(feRxSignal) || isempty(feRxLimits)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 5) || isempty(echoTimes)
    echoTimes = [];
end
if (nargin < 6) || isempty(reversePolarity)
    reversePolarity = 1;
end
if (nargin < 8) || isempty(feFullRwTime) || isempty(feFullRwSignal)
    feFullRwTime   = [];
    feFullRwSignal = [];
end
if (nargin < 9) || isempty(expControl)
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

% extract dimensions
lenFeSignal = numel(feTime); % encoding
lenFullFeRw = numel(feFullRwSignal); % full rewinder
% check correctness
if ( lenFullFeRw < 1 ) && (~reversePolarity)
    % ERROR: no way of reversing
    msg = sprintf( ' There is no RW to apply.' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
% number of time points and echo times and delays
numEchoes   = numel(echoTimes);
if reversePolarity
    numTimePoints   = (numEchoes+1)*lenFeSignal;
    minEchoTime     = feTime(end);
else
    % also apply RW
    numTimePoints   = (numEchoes+1)*lenFeSignal + numEchoes*lenFullFeRw;
    minEchoTime     = feTime(end) + feFullRwTime(end);
end
% inter echo delays
interEchoDelay = echoTimes - minEchoTime;
if interEchoDelay < 0
    % ERROR: not echoTimes too short
    msg = sprintf( ' Required Echo times too short, minimum %.3f ms.', minEchoTime*1e3 );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
% add time points in case we need interEchoDelay
numTimePoints = numTimePoints + nnz(interEchoDelay > 0);

% number of time points:
signalFE    = zeros(numTimePoints,1);
signalRX    = zeros(numTimePoints,1);
time        = zeros(numTimePoints,1);

% first FE encoding
targetIdx           = 1:lenFeSignal;
signalFE(targetIdx)	= feSignal(:);
signalRX(targetIdx) = feRxSignal(:);
time(targetIdx)     = feTime;
rxLimits            = feRxLimits; % shift the rxLimits indexes

% rest of the echoes
idxShift    = lenFeSignal;    % shift to last position
timeShift   = time(idxShift); % time shift
fePolarity  = 1;
for ii = 1:numEchoes
    % if same polarity, add full FE RW
    if reversePolarity
        fePolarity = -1*fePolarity;
    else
        % add full rewind
        targetIdx = (1:lenFullFeRw) + idxShift;
        signalFE(targetIdx) = -feFullRwSignal(:);
        time(targetIdx)     =  feFullRwTime(:) + timeShift;
        % update shifts
        idxShift    = idxShift + lenFullFeRw;
        timeShift   = time(idxShift);
    end
    % add any delay to correct time between echoes
    if interEchoDelay(ii) > 0
        timeShift           = timeShift + interEchoDelay(ii);
        idxShift            = idxShift + 1;
        time(idxShift)      = timeShift;
        signalFE(idxShift)  = 0;
    end
    % add new FE encoding
    signalFE( (1:lenFeSignal) + idxShift)   =  fePolarity*feSignal(:);
    signalRX( (1:lenFeSignal) + idxShift)   =  feRxSignal(:);
    time( (1:lenFeSignal) + idxShift)       =  feTime + timeShift;
    % append rxLimits, with correct shift
    rxLimits    = [rxLimits; feRxLimits + idxShift];
    % update shifts
    idxShift    = idxShift + lenFeSignal;
    timeShift   = time(idxShift);
end

if idxShift < numTimePoints
    % ERROR
    msg = sprintf( ' Mismatch in number of time points.');
    ME = MException(['warning:',functionName], '%s', msg);
    throw(ME);
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


