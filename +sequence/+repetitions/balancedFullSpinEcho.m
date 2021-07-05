function [repetition] = balancedFullSpinEcho(...
    acqData, mainRF, refRF, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.BALANCEDFULLSPINECHO
%
%	Generates Spin Echo repetition.
%
% INPUT
%   acqData
%   mainRF
%   refRF
%   mrSystem        
%   expControl      
%
% OUTPUT
%   repetition   repetition struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:repetitions:balancedFullSpinEcho';

% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFlie,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% generate refocusing (180+Encoding) Spin Echo repetition encoding
[repetition180] = sequence.repetitions.refocusSpinEcho(...
    acqData, refRF, mrSystem, expControl);

% if user forced min TE
if acqData.forceMinTE
    userTE = repetition180.TE;
else
    % user provided TE
    userTE  = acqData.TE;
end

%% generate the preparation block (90+FERW)

% update TE to be that provided by refocusing rep
acqData.TE = repetition180.TE;
effTE = repetition180.TE;

% pre-encoding gradient
areaFE = sum(repetition180.feSignal.*...
    [repetition180.time(1); diff(repetition180.time)]);

% generate the RF 90 Prep repetition 
[repetition90] = sequence.repetitions.prepSpinEcho(areaFE,...
    acqData, mainRF, mrSystem, expControl);

%% get slice selection gradient level in the 180 
ssGradLevel = repetition180.ssGradLevel;

%% verify feasibility and correctness of timing
% TE 
minTE = max(repetition180.TE, repetition90.TE);
if userTE < minTE
    msg = sprintf( ['The selected Echo Time (TE=%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is TE=%.2fms'],...
        userTE*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end

% TR
minTR    = repetition90.time(end) + repetition180.time(end);
if acqData.forceMinTR
    % force to minimum TR possible
    tWaitRep = 0;
    acqData.TR = minTR;
else
    tWaitRep = acqData.TR - minTR;
end
if tWaitRep < 0
    msg = sprintf( ['The selected Repetition Time (TR=%.2fms) is too short. ',...
        'Minimum Repetition Time for the current configuration is TR=%.2fms'],...
        acqData.TR*1e3, minTR*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end

%% combine the two repetitions
numSteps = repetition90.numSteps + repetition180.numSteps;
idx90  = 1:repetition90.numSteps;
idx180 = repetition90.numSteps+1:numSteps;

%% allocate space
repetition.time      = zeros(numSteps,1);
repetition.timeDiff  = zeros(numSteps,1);
repetition.rxSignal  = zeros(numSteps,1); % receiver readout
repetition.swcSignal = zeros(numSteps,1); % software crusher
repetition.feSignal  = zeros(numSteps,1); % freq encoding (x) gradient
repetition.peSignal  = zeros(numSteps,1); % phase encoding (y) gradient
repetition.seSignal  = zeros(numSteps,1); % 3D phase (slice) encoding (z) gradient
repetition.rfmSignal = zeros(numSteps,1); % RF magnitude
repetition.rfpSignal = zeros(numSteps,1); % RF phase
repetition.rffSignal = zeros(numSteps,1); % RF frequency
repetition.ssSignal  = zeros(numSteps,1); % slice selection gradient
repetition.rfEntries = zeros(numSteps,1); % ones or zeros depending on RF

%% assign blocks
repetition.time(idx90)  = repetition90.time(:);
repetition.time(idx180) = repetition180.time(:) + repetition90.totalTime;
repetition.time(end)    = acqData.TR;
repetition.timeDiff(:)  = diff([0; repetition.time]);
%% Pre-Encoding
repetition.rfEntries(idx90)	= repetition90.rfEntries;
repetition.rfmSignal(idx90)	= repetition90.rfmSignal;
repetition.rfpSignal(idx90)	= repetition90.rfpSignal;
repetition.ssSignal(idx90)	= repetition90.ssSignal;
repetition.feSignal(idx90)  = repetition90.feSignal;
repetition.peSignal(idx90)  = repetition90.peSignal;
repetition.seSignal(idx90)  = repetition90.seSignal;
repetition.rxSignal(idx90)  = repetition90.rxSignal;
%% Refocusing
repetition.rfEntries(idx180)    = repetition180.rfEntries;
repetition.rfmSignal(idx180)	= repetition180.rfmSignal;
repetition.rfpSignal(idx180)	= repetition180.rfpSignal;
repetition.ssSignal(idx180)     = repetition180.ssSignal;
repetition.feSignal(idx180)     = repetition180.feSignal;
repetition.peSignal(idx180)     = repetition180.peSignal;
repetition.seSignal(idx180)     = repetition180.seSignal;
repetition.rxSignal(idx180)     = repetition180.rxSignal;

%% pulse sequences info
repetition.type         = 'Spin-Echo';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = ssGradLevel;

%% for splitting into parts (RF vs Non-RF)
repetition.numParts     = repetition90.numParts + repetition180.numParts;
repetition.partLimits   = [repetition90.partLimits;...
    repetition180.partLimits + repetition90.numSteps];
repetition.partType     = repetition90.partType;
repetition.partType(repetition90.numParts+1:repetition.numParts) = ...
    repetition180.partType;

%% readout limits
repetition.rxLimits   = repetition180.rxLimits + repetition90.numSteps;

% update TE / TR
repetition.TE         = acqData.TE;
repetition.TR         = repetition.time(end);
% add info of minimum TR (including waiting for correct TE)
repetition.minTR      = minTR;

%% verify echo times
idxRF90     = find(repetition90.rfEntries>0);
idxRF180    = find(repetition180.rfEntries>0);
center90   = mean(repetition.time(idxRF90));
center180  = mean(repetition.time(idxRF180 + repetition90.numSteps));
T90to180   = center180 - center90;
T180toEcho = mean(repetition.time(repetition.rxLimits))- center180;
if abs(T90to180 - repetition.TE/2 ) > 1e-12
    ME = MException('sequence:wrongEchoTime',...
        '%s : TE and time between RFs do not match',functionName);
    throw(ME);
end
if abs(T180toEcho - repetition.TE/2 ) > 1e-12
    ME = MException('sequence:wrongEchoTime',...
        '%s : TE missaligned with center of readout',functionName);
    throw(ME);
end

%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s repetition, elapsed time %.3fs',...
        functionName, repetition.type, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', repetition.totalTime);
    fprintf(fid, '\n  Number steps        %d', repetition.numSteps);
    fprintf(fid, '\n  Number readouts     %d', repetition.numRxs);
    fprintf(fid, '\n  Effective TE        %.3fms', 1e3*repetition.TE );
    fprintf(fid, '\n  Effective TR        %.3fms', 1e3*repetition.TR );
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
%         % plot
%         figure();
%         plot(repetition.time,repetition.rfmSignal*1e3)
%         hold on
%         plot(repetition.time(rfLimits90),[0,0], '^');
%         plot(repetition.time(rfLimits180),[0,0], 'v');
%         plot(repetition.time,repetition.feSignal);
%         plot(repetition.time,repetition.peSignal);
%         plot(repetition.time,repetition.ssSignal);
%         plot(repetition.time(repetition.rxSignal>0), ...
%             repetition.feSignal(repetition.rxSignal>0), 'o')

end


