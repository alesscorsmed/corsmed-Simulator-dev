function [repetition] = fullSpinEcho(...
    acqData, mainRF, refRF, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.FULLSPINECHO
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
functionName = 'sequence.repetitions.fullSpinEcho';

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

%% generate unbalanced SE encoding
[encTime,signalFE,signalPE,signalSE,signalRX,areaFE,dtEnc,effTE,rxLimits] = ...
    sequence.encodings.cartesianPhaseBalanced(...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
    acqData.numFE, acqData.numPE, acqData.numSE,...
    acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
    expControl.sequence.dtGR, expControl );

%% generate Pre-Encoding Gradient
[rwTime,rwSignal,~,~,~] = sequence.waveforms.grTrapArea(...
        areaFE/2,mrSystem.maxGStrenght,mrSystem.SlewRate,dtEnc);

%% generate the waveforms for the RF 90 (use main RF info)
[rf90Time,rfm90,rfp90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% generate the waveforms for the RF 180 (use refocusing RF info)
[rf180Time,rfm180,rfp180,grSlice180,rfLimits180] = ...
    sequence.excitations.sliceSelectRF( refRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% get slice selection gradient level in the 180 
ssGradLevel = grSlice180(round(mean(rfLimits180)));

%% verify feasibility of TE/TR and correctness of timing
% TE
minTE90  = rf90Time(end)/2 + rwTime(end) + rf180Time(end)/2;
minTE180 = rf180Time(end)/2 + effTE;
minTE = 2*max(minTE180,minTE90);
if acqData.forceMinTE
    % force to minimum TE possible
    acqData.TE = minTE;
end
% for spin echo sequences, TE is from center of 90 to center of readout
TE = acqData.TE/2; % TE is from center of 180 to center of readout
tWait90  = TE - minTE90;
tWait180 = TE - minTE180;
% TR
minTR    = rf90Time(end)/2 + TE + rf180Time(end)/2 + tWait180 + encTime(end);
if acqData.forceMinTR
    % force to minimum TR possible
    tWaitRep = 0;
    acqData.TR = minTR;
else
    tWaitRep = acqData.TR - minTR;
end
if tWait90 < 0
    if ~isempty(expControl.connLocalDB) ...
            && strcmpi(expControl.application,'edutool')
        msg = sprintf( ['The selected Echo Time (TE=%.2fms) is too short. ',...
            'Minimum Echo Time for the current configuration is TE=%.2fms'],...
            acqData.TE*1e3, minTE*1e3 );
        eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error', expControl.experimentID, expControl.pulseqID);
    else
        ME = MException('sequence:wrongEchoTime',...
            '%s : %.1fms TE too short (min 2x%.1fms for pre-encoding between RFs)',...
            functionName, 1e3*acqData.TE, 1e3*minTE90);
        throw(ME);
    end
else
    tWait90 = max(tWait90,1e-12);
end
if tWait180 < 0
    if ~isempty(expControl.connLocalDB) ...
            && strcmpi(expControl.application,'edutool')
        msg = sprintf( ['The selected Echo Time (TE=%.2fms) is too short. ',...
            'Minimum Echo Time for the current configuration is TE=%.2fms'],...
            acqData.TE*1e3, minTE*1e3 );
        eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error', expControl.experimentID, expControl.pulseqID);
    else
        ME = MException('sequence:wrongEchoTime',...
            '%s : %.1fms TE too short (min 2x%.1fms for alligned echo)',...
            functionName, 1e3*acqData.TE, 1e3*minTE180);
        throw(ME);
    end
else
    tWait180 = max(tWait180,1e-12);
end
if tWaitRep < 0
    if ~isempty(expControl.connLocalDB) ...
            && strcmpi(expControl.application,'edutool')
        msg = sprintf( ['The selected Repetition Time (TR=%.2fms) is too short. ',...
            'Minimum Repetition Time for the current configuration is TR=%.2fms'],...
            acqData.TR*1e3, minTR*1e3 );
        eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error', expControl.experimentID, expControl.pulseqID);
    else
        ME = MException('sequence:wrongRepTime',...
            '%s : %.1fms TR too short (minimum %.1fms)',...
            functionName, acqData.TR*1e3, minTR*1e3);
        throw(ME);
    end
else
    tWaitRep = max(tWaitRep,1e-12);
end

%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRF180   = length(rf180Time);
numTimesRw      = length(rwTime);
numTimesEnc     = length(encTime);

% start/end of RF90
iniRF90         = 1;
endRF90         = iniRF90 + numTimesRF90-1;
rfLimits90      = rfLimits90 + iniRF90-1;
% start/end of 1st gradient
idx1Wait90      = endRF90 + 1;
iniRw           = idx1Wait90 + 1;
endRw           = iniRw+numTimesRw-1;
idx2Wait90      = endRw + 1;
% start/end of RF180
iniRF180        = idx2Wait90 + 1;
endRF180        = iniRF180 + numTimesRF180-1;
rfLimits180     = rfLimits180 + iniRF180-1;
% start/end of Encoding gradients
idxWait180      = endRF180 + 1;
iniEnc          = idxWait180 + 1;
endEnc          = iniEnc+numTimesEnc-1;
idxWaitRep      = endEnc + 1;
% correct rxLimits
rxLimits        = rxLimits + idxWait180;
% total time of repetition
numSteps        = idxWaitRep;

%% pulse sequence

% waveforms
repetition.time      = zeros(numSteps,1);
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

% for area compressed sequences
repetition.hArea     = zeros(numSteps,1); % time deltas for areas
repetition.gxArea    = zeros(numSteps,1); % x gradient area delta
repetition.gyArea    = zeros(numSteps,1); % y gradient area delta
repetition.gzArea    = zeros(numSteps,1); % z gradient area delta

%% assign blocks
repetition.time(iniRF90:endRF90)     = rf90Time;
repetition.time(idx1Wait90)          = tWait90/2 + repetition.time(idx1Wait90-1);
repetition.time(iniRw:endRw)         = rwTime + repetition.time(iniRw-1);
repetition.time(idx2Wait90)          = tWait90/2 + repetition.time(idx2Wait90-1);
repetition.time(iniRF180:endRF180)   = rf180Time + repetition.time(iniRF180-1);
repetition.time(idxWait180)          = tWait180 + repetition.time(idxWait180-1);
repetition.time(iniEnc:endEnc)       = encTime + repetition.time(iniEnc-1);
repetition.time(idxWaitRep)          = tWaitRep + repetition.time(idxWaitRep-1);
repetition.timeDiff(:)               = diff([0; repetition.time]);
%% RF
repetition.rfEntries(rfLimits90(1):rfLimits90(2))   = 1;
repetition.rfmSignal(iniRF90:endRF90)    = rfm90;
repetition.rfpSignal(iniRF90:endRF90)    = rfp90;
repetition.ssSignal(iniRF90:endRF90)     = grSlice90;
%% Refocusing
repetition.rfEntries(rfLimits180(1):rfLimits180(2))   = 1;
repetition.rfmSignal(iniRF180:endRF180)  = rfm180;
repetition.rfpSignal(iniRF180:endRF180)  = rfp180;
repetition.ssSignal(iniRF180:endRF180)   = grSlice180;
%% Pre-Rewind
repetition.feSignal(iniRw:endRw) = ...
    repetition.feSignal(iniRw:endRw) + rwSignal(:);
%% Encoding
repetition.feSignal(iniEnc:endEnc) = ...
    repetition.feSignal(iniEnc:endEnc) + signalFE(:);
repetition.peSignal(iniEnc:endEnc) = ...
    repetition.peSignal(iniEnc:endEnc) + signalPE(:);
repetition.seSignal(iniEnc:endEnc) = ...
    repetition.seSignal(iniEnc:endEnc) + signalSE(:);
repetition.rxSignal(iniEnc:endEnc) = signalRX;

%% pulse sequences info
repetition.type         = 'Spin-Echo';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = ssGradLevel;

%% for splitting into parts (RF vs Non-RF)
if expControl.sequence.minContextExc
    % reduce context change in Kernel exec
    % by creating a unique part
    % caution: this may cause issues with numerical simulations
    numParts = 1;
    partLimits = [1, rfLimits180(2)];
    partType{numParts} = 'RF';
else
    % multiple parts, depending on its characteristics
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
end
% rest may have readouts
numParts = numParts + 1;
partLimits = [partLimits; rfLimits180(2)+1, repetition.numSteps];
partType{numParts} = 'RO';
% assign
repetition.numParts   = numParts; % number of parts
repetition.partType   = partType; % type of part
repetition.partLimits = partLimits; % index start/end of part
repetition.rxLimits   = rxLimits; % index start/end of readouts
% update TE / TR
repetition.TE         = acqData.TE;
repetition.TR         = repetition.time(end);

%% verify echo times
T90to180   = repetition.time(rfLimits180) - repetition.time(rfLimits90);
T180toEcho = mean(repetition.time(repetition.rxSignal>0)) ...
    - mean(repetition.time(rfLimits180));
if (T90to180 - repetition.TE/2 ) > 1e-12
    ME = MException('sequence:wrongEchoTime',...
        '%s : TE and time between RFs do not match',functionName);
    throw(ME);
end
if (T180toEcho - repetition.TE/2 ) > 1e-12
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


