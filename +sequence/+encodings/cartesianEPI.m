function [time,signalFE,signalPE,signalRX,dt,effTE,effESP,rxLimits,numREP] = ...
    cartesianEPI( receiverBW, fovFE, fovPE, numFE, numPE, ...
    samplingFactorFE, samplingFactorPE, ...
    readPartFourier, phasePartFourier, rFactor,...
    maxGStrenght, gradSlewrate, gamma, dt, expControl)
%
% SEQUENCE.ENCODINGS.CARTESIANEPI
%
%	Generates a cartesian EPI encoding.
%
% INPUT
%   receiverBW      bandwidth of HW receiver in Hz (will determine dt)
%   fovFE           Field of View in the frequency encoding direction
%   fovPE           Field of View in the phase encoding direction
%   numFE           number of frequency encodings (samples in readout)
%   numPE           number of phase encodings
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
%   signalRX        readout signal
%   dt              time step used
%   effTE           effective Time Echo: time to the center of the encoding
%   effESP          echo time spacing, time between readouts
%   rxLimits        start and end indexes of the readouts (numPE x 2)
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:encodings:cartesianEPI';

if (nargin < 1 || isempty(receiverBW))
    receiverBW = 200e3;  % in Hz (this is a feature of the hardware)
end
if (nargin < 2 || isempty(fovFE))
    fovFE   = 0.3;
end
if (nargin < 3 || isempty(fovPE))
    fovPE 	= 0.3;
end
if (nargin < 4 || isempty(numFE))
    numFE	= 256;
end
if (nargin < 5 || isempty(numPE))
    numPE	= 128;
end
if (nargin < 6 || isempty(samplingFactorFE))
    samplingFactorFE = 1;
end
if (nargin < 7 || isempty(samplingFactorPE))
    samplingFactorPE = 1;
end
if (nargin < 8 || isempty(readPartFourier))
    readPartFourier = 1.0;
end
if (nargin < 9 || isempty(phasePartFourier))
    phasePartFourier = 1.0;
end
if (nargin < 10 || isempty(rFactor))
    rFactor = 1;
end
if (nargin < 11 || isempty(maxGStrenght))
    maxGStrenght = 0.030;
end
if (nargin < 12 || isempty(gradSlewrate))
    gradSlewrate = 150.0;
end
if (nargin < 13 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end
if (nargin < 14 || isempty(dt))
    dt = 1e-6; % max allowed dt
end
if (nargin < 15 || isempty(expControl))
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
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end
% generate waveform
gradTime = numFE/simRxBW;
[feTime,feSignal,feArea,~,fePlatLimits] = ...
    sequence.waveforms.grTrapPlateau(gradTime,gradAmp,gradSlewrate,dt);

%% generate PE blips
blipArea = rFactor/(gamma*fovPE);
[peBlipTime,peBlipSignal,peBlipArea,~,~] = ...
    sequence.waveforms.grTrapArea(blipArea,maxGStrenght,gradSlewrate,dt);

%% construct single part of the EPI encoding (Freq Encoding + Blip)
% try to fit the blip during the ramps
lenBlip     = length(peBlipTime);
lenFeRise   = fePlatLimits(1);
lenFeSignal = length(feTime);
lenFeDelay  = lenBlip - lenFeRise + 1; % time we need to add to the FE
lenEnc      = lenFeSignal;
if lenFeDelay > 0
    % add extra time steps
    lenEnc = lenEnc + lenFeDelay; 
end
% assign encoding signals
feEnc  = zeros(lenEnc,1);
peEnc  = zeros(lenEnc,1);
roEnc  = zeros(lenEnc,1);
% assign signals
feEnc(1:lenFeSignal)            = feSignal(:);
peEnc(lenEnc-lenBlip+1:lenEnc)  = peBlipSignal(:);
% readout placement
roIdx = fePlatLimits(1):samplingRate:fePlatLimits(1)+samplingRate*numFE;
roIdx = roIdx(1:numFE);
roIdx = roIdx + floor( (fePlatLimits(2) - roIdx(end))/2 );
roEnc(roIdx) = 1:numFE;

%% repeat the encoding with the encoding scaling
% number of repetition depends on the partial Fourier
numREP = round(phasePartFourier*numPE);
% encoding levels
peEncodingLevels = ones(numREP,1);
feEncodingLevels = ones(numREP,1);
roEncodingLevels = ones(numREP,1);
% Parallel acceleration
peEncodingLevels = peEncodingLevels(1:rFactor:end,1);
feEncodingLevels = feEncodingLevels(1:rFactor:end,1);
roEncodingLevels = roEncodingLevels(1:rFactor:end,1);
% final number of encodings after acceleration
numENC = length(peEncodingLevels);
% EPI related
peEncodingLevels(end)       = 0; % do not apply last blip
feEncodingLevels(2:2:end)   = -1; % negative for even encodings
% apply encodings
feEnc = kron(feEncodingLevels,feEnc);
peEnc = kron(peEncodingLevels,peEnc);
roEnc = kron(roEncodingLevels,roEnc);

% index limits of the readouts
rxLimitsCol = lenEnc*(0:numENC-1).';
rxLimits    = [rxLimitsCol + roIdx(1), rxLimitsCol + roIdx(end)];

% Echo time spacing: duration of individual encoding
effESP  = lenEnc;
% Effective Echo Time
if phasePartFourier < 1.0
    preRwScale = (phasePartFourier - 0.5)*numPE;
    effTE = lenEnc*(numREP - ceil(numPE/2));
else
    preRwScale = 0.5*numPE;
    effTE = lenEnc*ceil(numENC/2);
end


%% generate FE gradient signal for Rewind
[feRwTime,feRwSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(feArea/2,0.95*maxGStrenght,gradSlewrate,dt);

%% generate PE gradient signal for Pre-Rewind
% notice that we rewind for partial Fourier Encoding
pePreRwArea  = preRwScale/(gamma*fovPE);
[pePreRwTime,pePreRwSignal,pePreRwArea,~,~] = ...
    sequence.waveforms.grTrapArea(pePreRwArea,0.95*maxGStrenght,gradSlewrate,dt);
% notice that we post rewind for all accrued area
pePostRwArea  = (numENC-1)*peBlipArea -pePreRwArea ;
[pePostRwTime,pePostRwSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(pePostRwArea,0.95*maxGStrenght,gradSlewrate,dt);

%% combine parts into signals
% last active entry of the frequency encoding
lenEnc          = find(abs(feEnc) > 0);
lenEnc          = lenEnc(end);
lenFeRw         = length(feRwTime);
lenPrePeRw      = length(pePreRwTime);
lenPostPeRw     = length(pePostRwTime);
lenPreRW        = max(lenFeRw,lenPrePeRw);
lenPostRW       = max(lenFeRw,lenPostPeRw);
numTimePoints   = lenEnc + lenPreRW + lenPostRW;
effTE           = effTE + lenPreRW;
rxLimits        = rxLimits + lenPreRW;
% create signals
signalFE = zeros(numTimePoints,1);
signalPE = zeros(numTimePoints,1);
signalRX = zeros(numTimePoints,1);
time     = zeros(numTimePoints,1);
% populate
signalFE(1:lenFeRw)                                     = -1*feRwSignal(:);
signalPE(1:lenPrePeRw)                                  = -1*pePreRwSignal(:);
signalFE(lenPreRW+1:lenPreRW+lenEnc)                    = feEnc(1:lenEnc);
signalPE(lenPreRW+1:lenPreRW+lenEnc)                    = peEnc(1:lenEnc);
signalRX(lenPreRW+1:lenPreRW+lenEnc)                    = roEnc(1:lenEnc);
signalFE(lenPreRW+lenEnc+1:lenPreRW+lenEnc+lenFeRw)     = (-1)^(numENC)*feRwSignal(:);
signalPE(lenPreRW+lenEnc+1:lenPreRW+lenEnc+lenPostPeRw) = -1*pePostRwSignal(:);
time(1:end)                                             = dt*(1:numTimePoints);

% apply partial echo if needed
partNumFE = round(readPartFourier*numFE);
startFE   = 1 + numFE - partNumFE;
% find the indexes of the positions of the Sampling Points
idxRx = find( signalRX > 0);
% reshape the indexes: numFE x numEchoes
idxRx = reshape(idxRx,[],numENC);
% get the new indexes by cropping the first startFE-1
newIdxRx = idxRx(startFE:end,:);
% zero all previous rx Signals
signalRX(:) = 0;
% assign new values to the partial acquisition entries
signalRX(newIdxRx(:)) = 1:length(newIdxRx(:));

% scale timing
effESP  = effESP*dt;
effTE   = effTE*dt;

% verify areas null out
if abs(sum(signalFE)*dt) > 1e-12    
    ME = MException('sequence:wrongGradArea',...
        '%s : FE gradient area is not zero after balanced encoding',functionName);
    throw(ME);
end
if abs(sum(signalPE)*dt) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : PE gradient area is not zero after balanced encoding',functionName);
    throw(ME);
end

% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', time(end));
    fprintf(fid, '\n  Echo Time Spacing   %.3fms', effESP*1e3 );
    fprintf(fid, '\n  Eff. Time Echo      %.3fms', effTE*1e3 );
    fprintf(fid, '\n  Rx Bandwidth        %.3fKHz', simRxBW*1e-3);
    fprintf(fid, '\n  Rx Sampling Rate    %.3fus', 1e6/simRxBW);
    fprintf(fid, '\n  Time step           %.3fus', dt*1e6);
    fprintf(fid, '\n  Number steps        %d', numTimePoints);
    fprintf(fid, '\n  Number readouts     %d (%d per encoding)', nnz(signalRX), nnz(signalRX)/numENC);
    fprintf(fid, '\n  Number encodings    %d', numENC );
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
    %     % plot
    %     plot(time, signalFE)
    %     hold on
    %     plot(time, signalPE)
    %     hold on
    %     plot(time(signalRX>0), signalFE(signalRX>0), 'o')
    
end

