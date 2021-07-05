function [rwTime,rwSignalFE,rwSignalPE,rwSignalSE,...
    roTime, roSignalFE, roSignalRX, roAreaFE, dt, effTE, rxLimits] = ...
    cartesianGeneral( receiverBW, fovFE, fovPE, fovSE,...
    numFE, numPE, numSE, ...
    samplingFactorFE, samplingFactorPE, samplingFactorSE, ...
    maxGStrenght, gradSlewrate, ...
    gamma, dt, expControl)
%
% SEQUENCE.ENCODINGS.CARTESIANGENERAL
%
%	Generates a cartesian encoding with the encoding on one side:
%   FE:  HalfRewind
%   PE:  HalfRewind 
%   SE:  HalfRewind 
%   and the readout in the other
%   FE:  Encoding 
%   each on independent time 
%
% INPUT
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
functionName = 'sequence.encodings.cartesianPhaseBalanced';

if (nargin < 1 || isempty(receiverBW))
    receiverBW = 200e3;  % in Hz (this is a feature of the hardware)
end
if (nargin < 2 || isempty(fovFE))
    fovFE   = 0.3;
end
if (nargin < 3 || isempty(fovPE))
    fovPE 	= 0.3;
end
if (nargin < 4 || isempty(fovSE))
    fovSE 	= 0.3;
end
if (nargin < 5 || isempty(numFE))
    numFE	= 256;
end
if (nargin < 6 || isempty(numPE))
    numPE	= 128;
end
if (nargin < 7 || isempty(numSE))
    numSE	= 128;
end
if (nargin < 8 || isempty(samplingFactorFE))
    samplingFactorFE = 1;
end
if (nargin < 9 || isempty(samplingFactorPE))
    samplingFactorPE = 1;
end
if (nargin < 10 || isempty(samplingFactorSE))
    samplingFactorSE = 1;
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
    expControl.debug.debugMode = 0;
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
[roTime,roSignalFE,roAreaFE,~,roPlatLimits] = ...
    sequence.waveforms.grTrapPlateau(gradTime,gradAmp,gradSlewrate,dt);

%% generate FE rewinder
[feTime,feSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(roAreaFE/2,maxGStrenght,gradSlewrate,dt);

%% generate half of PE encoding
peArea = numPE/(gamma*fovPE)/2;
[peTime,peSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(peArea,maxGStrenght,gradSlewrate,dt);

%% generate half of 3D encoding
seArea = numSE/(gamma*fovSE)/2;
[seTime,seSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(seArea,maxGStrenght,gradSlewrate,dt);

%% construct single part of the RW encoding
lenFe       = length(feTime);
lenPe       = length(peTime);
lenSe       = length(seTime);
lenRW       = max([lenFe,lenPe,lenSe]);
% create signals
rwSignalFE = zeros(lenRW,1);
rwSignalPE = zeros(lenRW,1);
rwSignalSE = zeros(lenRW,1);
rwTime     = zeros(lenRW,1);
% populate
rwSignalFE(1:lenFe) =  feSignal(:);
rwSignalPE(1:lenPe) =  peSignal(:);
rwSignalSE(1:lenSe) =  seSignal(:);
rwTime(1:lenRW)     = dt*(1:lenRW);

% readout placement
roIdx = roPlatLimits(1):samplingRate:roPlatLimits(1)+samplingRate*numFE;
roIdx = roIdx(1:numFE);
roIdx = roIdx + floor( (roPlatLimits(2) - roIdx(end))/2 );
roSignalRX      = 0*roSignalFE;
roSignalRX(roIdx) = 1:numFE;
% index limits of first readout
rxLimits = [roIdx(1), roIdx(end)]; 
% effective echo center of RX
effTE  = mean(roTime(rxLimits));

% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  Eff. Time Echo      %.3fms', effTE*1e3 );
    fprintf(fid, '\n  Rx Bandwidth        %.3fKHz', simRxBW*1e-3);
    fprintf(fid, '\n  Rx Sampling Rate    %.3fus', 1e6/simRxBW);
    fprintf(fid, '\n  Time step           %.3fus', dt*1e6);
    fprintf(fid, '\n  Number readouts     %d', nnz(roSignalRX));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
    %     % plot
    %     plot(time, signalFE)
    %     hold on
    %     plot(time, signalPE)
    %     hold on
    %     plot(time, signalSE)
    %     hold on
    %     plot(time(signalRX>0), signalFE(signalRX>0), 'o')
    
end

