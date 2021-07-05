function [time, signalFE, signalPE, signalSE,...
    signalRX, rxLimits, effTE, feArea, dt] = ...
    cartesian( receiverBW, ...
    fovFE, fovPE, fovSE, numFE, numPE, numSE, ...
    samplingFactorFE, samplingFactorPE, samplingFactorSE, ...
    freqBalance, phaseBalance, sliceBalance, ...
    echoTimes, reverseEchoPolarity,... 
    maxGStrenght, gradSlewrate, gamma, dt, expControl)
%
% SEQUENCE.ENCODINGS.CARTESIAN
%
%	Generates a general cartesian encoding.
%   The *Balance parameters define if we apply pre and/or post rewind.
%   The echoTimes and Polarity define the multi-Echo pattern
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
functionName = 'sequence:encodings:cartesian';
if (nargin < 1 || isempty(receiverBW))
    receiverBW = 200e3;  % in Hz (this is a feature of the hardware)
end
if (nargin < 2) || isempty(fovFE)
    fovFE   = 0.3;
end
if (nargin < 3) || isempty(fovPE)
    fovPE 	= 0.3;
end
if (nargin < 4) || isempty(fovSE)
    fovSE 	= 0.3;
end
if (nargin < 5) || isempty(numFE)
    numFE	= 256;
end
if (nargin < 6) || isempty(numPE)
    numPE	= 128;
end
if (nargin < 7) || isempty(numSE)
    numPE	= 128;
end
if (nargin < 8) || isempty(samplingFactorFE)
    samplingFactorFE = 1;
end
if (nargin < 9) || isempty(samplingFactorPE)
    samplingFactorPE = 1;
end
if (nargin < 10) || isempty(samplingFactorSE)
    samplingFactorSE = 1;
end
if (nargin < 11) || isempty(freqBalance)
    freqBalance = 'full';
end
if (nargin < 12) || isempty(phaseBalance)
    phaseBalance = 'full';
end
if (nargin < 13) || isempty(sliceBalance)
    sliceBalance = 'full';
end
if (nargin < 14) || isempty(echoTimes)
    echoTimes = []; % no echo
end
if (nargin < 15) || isempty(reverseEchoPolarity)
    reverseEchoPolarity = 1;
end
if (nargin < 16) || isempty(maxGStrenght)
    maxGStrenght = 0.030;
end
if (nargin < 17) || isempty(gradSlewrate)
    gradSlewrate = 150.0;
end
if (nargin < 18) || isempty(gamma)
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end
if (nargin < 19) || isempty(dt)
    dt = 1e-6; % max allowed dt
end
if (nargin < 20) || isempty(expControl)
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
dtBW    = 1/simRxBW;       % sampling rate during readout
dt      = min(dtBW,dt);    % time resolution

%% generate FE gradient signal for Readout
% check if maximum gradient limit is trespassed in Readout
gradAmp  = simRxBW/(gamma*fovFE);
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

% modify the waveform to force the samples at sampling rate
[feTime, feSignal, fePlatLimits, rxSignal, rxLimits] = ...
    sequence.tools.readoutPlateauPlacement( ...
    numFE, dtBW, feTime, feSignal, fePlatLimits, expControl );

% get the number of steps for the fe Ramp
feRiseSteps = fePlatLimits(1);

%% Prepare echo info 
numEchoes = numel(echoTimes);
if numEchoes > 0
    % minimum echo time
    minEchoTime = feTime(end);
    % check if we need to mantain the echo polarity
    if ~reverseEchoPolarity
        % we need to Full RW each encoding
        [feFullRwTime,feFullRwSignal,~,~,~] = ...
            sequence.waveforms.grTrapArea(feArea,maxGStrenght,gradSlewrate,dt);
        minEchoTime = minEchoTime + feFullRwTime(end);
    else
        feFullRwTime    = [];
        feFullRwSignal  = [];
    end
    % verify we can fit the echo times
    if minEchoTime > min(echoTimes)
        msg = sprintf( ['Times between Echoes (%.3fms) '...
            'are too short (minimum %.3fms). '...
            'Consider increasing time between echoes.'],...
            min(echoTimes)*1e3, minEchoTime*1e3 );
        if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
            eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
        else
            ME = MException(['userwarning:',functionName], '%s', msg);
            throw(ME);
        end
    end
    % assemble a fe encoding with all echoes
    [feTime, feSignal, rxSignal, rxLimits] = ...
        sequence.tools.multiEchoAssembly( ...
        feTime, feSignal, rxSignal, rxLimits, ...
        echoTimes, reverseEchoPolarity, feFullRwTime, feFullRwSignal, ...
        expControl );    
end



SpoilAmpl = 1.1;


%% generate FE rewinder
%  check the balancing of the FE encoding
switch lower(freqBalance)
    case 'full'
        feRwScale = [-1, -1];
        % in case of reversed echos, + balancing for even number of echoes
        if reverseEchoPolarity
            feRwScale(2) = (-1)^(numEchoes+1);
        end
    case 'prep'
        feRwScale = [-1,  0];
    case 'post'
        feRwScale = [ 0, -1];
    case 'spoil'
        feRwScale = [-1, SpoilAmpl];
        % in case of reversed echos, + reverse for even number of echoes
        if reverseEchoPolarity
            feRwScale(2) = -SpoilAmpl*(-1)^(numEchoes+1);
        end
    otherwise
        feRwScale = [ 0,  0];
end
% if needed, generate signal
if nnz(feRwScale)
    [~,feRwSignal,~,~,~] = ...
        sequence.waveforms.grTrapArea(feArea/2,maxGStrenght,gradSlewrate,dt);
else
    feRwSignal = [];
end

%% generate half of PE encoding
%  check the balancing of the PE encoding
switch lower(phaseBalance)
    case 'full'
        peRwScale = [ 1, -1];
    case 'prep'
        peRwScale = [ 1,  0];
    case 'post'
        peRwScale = [ 0, -1];
    case 'spoil'
        peRwScale = [ 1, SpoilAmpl];
    otherwise
        peRwScale = [ 0,  0];
end
% if needed, generate signal
peArea = numPE/(gamma*fovPE)/2;
if nnz(peRwScale)
    [~,peRwSignal,~,~,~] = ...
        sequence.waveforms.grTrapArea(peArea,maxGStrenght,gradSlewrate,dt);
else
    peRwSignal = [];
end

%% generate half of 3D encoding
%  check the balancing of the SE encoding
switch lower(sliceBalance)
    case 'full'
        seRwScale = [ 1, -1];
    case 'prep'
        seRwScale = [ 1,  0];
    case 'post'
        seRwScale = [ 0, -1];
    otherwise
        seRwScale = [ 0,  0];
end
% if needed, generate signal
seArea = numSE/(gamma*fovSE)/2;
if nnz(seRwScale)
    [~,seRwSignal,~,~,~] = ...
        sequence.waveforms.grTrapArea(seArea,maxGStrenght,gradSlewrate,dt);
else
    seRwSignal = [];
end

%% assemble the encoding signals with the RW
[time, signalFE, signalPE, signalSE, signalRX, rxLimits] = ...
    sequence.tools.encodingAssembly( feTime, feSignal, rxSignal, rxLimits, ...
    feRiseSteps, feRwSignal, peRwSignal, seRwSignal, ...
    feRwScale, peRwScale, seRwScale, dt, expControl );

% compute effective Echo time and assign RX indexes
effTE = mean(time(rxLimits(1,:)));
signalRX(signalRX>0) = 1:nnz(signalRX);
timeDiff = reshape(diff([0; time]),[],1);

%% verify areas null out depending on the balancing
% Freq
totalAreaFE = feArea + sum(feRwScale)*feArea/2;
if (numEchoes > 0) && reverseEchoPolarity
    totalAreaFE = totalAreaFE -rem(numEchoes,2)*feArea;
end
if abs(sum(timeDiff.*signalFE) - totalAreaFE) > 1e-12
    msg = sprintf( 'FE gradient total area mismatch with required balance scheme.' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
% Phase
totalAreaPE = sum(peRwScale)*peArea;
if abs(sum(timeDiff.*signalPE) - totalAreaPE) > 1e-12
    msg = sprintf( 'PE gradient total area mismatch with required balance scheme.' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end
% Slice
totalAreaSE = sum(seRwScale)*seArea;
if abs(sum(timeDiff.*signalSE) - totalAreaSE) > 1e-12
    msg = sprintf( 'SE gradient toal area mismatch with required balance scheme.' );
    ME = MException(['error:',functionName], '%s', msg);
    throw(ME);
end

%% check alignment of echo time and zero FE area
if strcmpi(freqBalance, 'full') || strcmpi(freqBalance, 'prep')
    echoAreaFE = 0;
else
    echoAreaFE = feArea/2;
end
if rem(numFE,2) % odd number of encodings: sample at echo
    idxEcho = rxLimits(1) + floor((numFE-1)/2);
    if abs(sum(timeDiff(1:idxEcho).*signalFE(1:idxEcho)) - echoAreaFE) > 1e-12
        msg = sprintf( 'Sampling placement mismatch at echo time.' );
        ME = MException(['error:',functionName], '%s', msg);
        throw(ME);
    end
else
    % find indexes before and after echo
    idxEcho1 = rxLimits(1) + floor((numFE-1)/2);
    idxEcho2 = idxEcho1 + 1;
    areaEcho1 = sum(timeDiff(1:idxEcho1).*signalFE(1:idxEcho1));
    areaEcho2 = sum(timeDiff(1:idxEcho2).*signalFE(1:idxEcho2));
    if abs( (areaEcho1 + areaEcho2)/2 - echoAreaFE) > 1e-12 && lower(freqBalance) ~= "spoil"
        msg = sprintf( 'Sampling placement not symmetric around echo time.' );
        ME = MException(['error:',functionName], '%s', msg);
        throw(ME);
    end
end

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
    fprintf(fid, '\n  Number steps        %d', numel(time));
    fprintf(fid, '\n  Number readouts     %d', nnz(signalRX));
    fprintf(fid, '\n  Number Echoes       %d', numEchoes);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
%     figure();
%     
%     figAx(1) = subplot(4,1,1);
%     plot(time, signalFE); hold on
%     %plot(newTime, newSignalFE); hold on
%     
%     figAx(2) = subplot(4,1,2);
%     plot(time, signalPE); hold on
%     %plot(newTime, newSignalPE); hold on
%     
%     figAx(3) = subplot(4,1,3);
%     plot(time, signalSE); hold on
%     %plot(newTime, newSignalSE); hold on
%     
%     figAx(4) = subplot(4,1,4);
%     plot(time(signalRX>0), signalFE(signalRX>0),'o'); hold on
%     %plot(newTime(newSignalRX>0), newSignalFE(newSignalRX>0), '+'); hold on
%     
%     linkaxes(figAx,'x');
%     
end

