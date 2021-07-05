function [repetition] = prepSpinEcho(areaFE,...
    acqData, mainRF, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.PREPSPINECHO
%
%	Generates preparation (90+FERW) for Spin Echo sequences.
%
% INPUT
%   acqData   
%   mainRF
%   mrSystem        
%   expControl      
%
% OUTPUT
%   repetition   repetition struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.repetitions.prepSpinEcho';

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

%% generate Pre-Encoding Gradient
[rwTime,rwSignal,~,~,~] = sequence.waveforms.grTrapArea(...
        areaFE/2,mrSystem.maxGStrenght,mrSystem.SlewRate, ...
        expControl.sequence.dtGR );

%% generate the waveforms for the RF 90 (use main RF info)
[rf90Time,rfm90,rfp90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl);

%% get slice selection gradient level in the 180 
ssGradLevel = grSlice90(round(mean(rfLimits90)));

%% verify feasibility of TE/TR and correctness of timing
minTE90  = rf90Time(end) + rwTime(end);
tWait90  = acqData.TE/2 - minTE90;
if tWait90 < 0
    if ~isempty(expControl.connLocalDB) ...
            && strcmpi(expControl.application,'edutool')
        msg = sprintf( ['The selected Echo Time (TE=%.2fms) is too short. ',...
            'Minimum Echo Time for the current configuration is TE=%.2fms'],...
            acqData.TE*1e3, 2*minTE90*1e3 );
        eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error', expControl.experimentID, expControl.pulseqID);
    else
        ME = MException('sequence:wrongEchoTime',...
            '%s : %.1fms TE too short (min %.1fms for pre-encoding between RFs)',...
            functionName, 1e3*acqData.TE, 2*1e3*minTE90);
        throw(ME);
    end
else
    tWait90 = max(tWait90,1e-12);
end

%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRw      = length(rwTime);

% start/end of RF90
iniRF90         = 1;
endRF90         = iniRF90 + numTimesRF90-1;
rfLimits90      = rfLimits90 + iniRF90-1;
% start/end of 1st gradient
idx1Wait90      = endRF90 + 1;
iniRw           = idx1Wait90 + 1;
endRw           = iniRw+numTimesRw-1;
idx2Wait90      = endRw + 1;
% total time of repetition
numSteps        = idx2Wait90;

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
%% RF
repetition.rfEntries(rfLimits90(1):rfLimits90(2))   = 1;
repetition.rfmSignal(iniRF90:endRF90)    = rfm90;
repetition.rfpSignal(iniRF90:endRF90)    = rfp90;
repetition.ssSignal(iniRF90:endRF90)     = grSlice90;
%% Pre-Rewind
repetition.feSignal(iniRw:endRw) = ...
    repetition.feSignal(iniRw:endRw) + rwSignal(:);

%% pulse sequences info
repetition.type         = 'prep-Spin-Echo';
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
    partLimits = [1, repetition.numSteps];
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
    partLimits = [partLimits; rfLimits90(2)+1, repetition.numSteps];
    partType{numParts} = 'GR';
end
% assign
repetition.numParts   = numParts; % number of parts
repetition.partType   = partType; % type of part
repetition.partLimits = partLimits; % index start/end of part
repetition.rxLimits   = []; % index start/end of readouts
% update TE / TR
repetition.TE         = acqData.TE;
repetition.TR         = repetition.time(end);

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


