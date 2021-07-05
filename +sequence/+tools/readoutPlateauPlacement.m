function [newTime, newSignalFE, newPlatLimits, newSignalRX, newRxLimits] = ...
    readoutPlateauPlacement( numRx, dtBW, time, signalFE, platLimits, expControl )
%
% SEQUENCE.TOOLS.READOUTPLACEMENT
%
%	Applies readout sampling given by a sampling rate, 
%   centered in an interval defined by the plateau limits.
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
functionName = 'sequence.tools.readoutPlateauPlacement';

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

% echo time, and symmetric readouts around it
effTE       = time(end)/2;
roTimeIni   = effTE - (numRx-1)/2*dtBW;
roTimeEnd   = roTimeIni + (numRx-1)*dtBW;
roTime      = roTimeIni:dtBW:roTimeEnd;
% plateau time limits
%  notice we each sample is associated with a time period that ends at the
%  time point!
if platLimits(1) == 1
    plTimeIni = 0.0;
else
    plTimeIni = time(platLimits(1)-1);
end
plTimeEnd = time(platLimits(2));
% check that readout does not fall outside plateau
if (roTimeIni < plTimeIni ) || (roTimeEnd > plTimeEnd)
    ME = MException('sequence:wrongReadoutPlacement',...
        '%s : readout times (%e - %e) out of plateau range (%e - %e)',...
        functionName, roTimeIni, roTimeEnd, ...
        time(platLimits(1)), time(platLimits(2)) );
    throw(ME);
end
% indexes outside the plateau to keep
idxKeep = union(1:platLimits(1)-1,platLimits(2):numel(time));
oldTime = time(idxKeep);
% new time
newTime = reshape(union(oldTime, roTime),[],1);
numTime = numel(newTime);
% get true/false depending of whether time is RX or not
idxRx   = ismember(newTime,roTime);
idxOld  = ismember(newTime,oldTime);
% allocate new signal arrays
newSignalFE         = zeros(numTime,1);
newSignalRX         = zeros(numTime,1);
% assign signals
newSignalFE(idxRx)  = signalFE(platLimits(1));
newSignalFE(idxOld) = signalFE(idxKeep);
newSignalRX(idxRx)  = 1:numRx;
% new limits Rx
idxRx               = find(idxRx);
newRxLimits         = [min(idxRx), max(idxRx)];
% new plateau Limits on new Time frame
idxPlEnd = find(ismember(newTime,plTimeEnd)); % index of Plateau End
if platLimits(1) == 1
    idxPlIni = 1;
else
    idxPlIni = find(ismember(newTime,plTimeIni))+1; % index of Plateau Start
end
newPlatLimits  = [idxPlIni, idxPlEnd];


%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  # time points %d', numTime);
    fprintf(fid, '\n  # samples     %d', numRx);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
