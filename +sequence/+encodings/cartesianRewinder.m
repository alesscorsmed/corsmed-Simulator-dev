function [ rwTime,rwSignalFE,rwSignalPE,rwSignalSE ] = ...
    cartesianRewinder( areaFE, areaPE, areaSE, ...
    maxGStrenght, gradSlewrate, dt, expControl)
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
functionName = 'sequence:encodings:cartesianRewinder';

if (nargin < 1 || isempty(areaFE))
    areaFE   = 0.0;
end
if (nargin < 2 || isempty(areaPE))
    areaPE 	= 0.0;
end
if (nargin < 3 || isempty(areaSE))
    areaSE 	= 0.0;
end
if (nargin < 4 || isempty(maxGStrenght))
    maxGStrenght = 0.030;
end
if (nargin < 5 || isempty(gradSlewrate))
    gradSlewrate = 150.0;
end
if (nargin < 6 || isempty(dt))
    dt = 1e-6; % max allowed dt
end
if (nargin < 7 || isempty(expControl))
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

%% generate FE rewinder
[feTime,feSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(areaFE,maxGStrenght,gradSlewrate,dt);

%% generate half of PE encoding
[peTime,peSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(areaPE,maxGStrenght,gradSlewrate,dt);

%% generate half of 3D encoding
[seTime,seSignal,~,~,~] = ...
    sequence.waveforms.grTrapArea(areaSE,maxGStrenght,gradSlewrate,dt);

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


% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  RW duration     %.3fus', rwTime(end)*1e6);
    fprintf(fid, '\n  Number steps    %d', lenRW);
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

