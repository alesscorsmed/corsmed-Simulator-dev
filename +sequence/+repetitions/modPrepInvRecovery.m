function [repetition] = prepInvRecovery(...
    acqData, prepIR, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.PREPINVRECOVERY
%
%	Generates preparation (180+TI) for IR sequences.
%
% INPUT
%   acqData   
%   prepIR
%   mrSystem        
%   expControl      
%
% OUTPUT
%   repetition   repetition struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:repetitions:prepInvRecovery';

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

%% generate the waveforms for the RF 90 (use prepIR info)
[rfIRTime,rfmIR,rfpIR,rffIR,grSliceIR,rfLimitsIR] = ...
    sequence.excitations.sliceSelectRF( prepIR,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% get slice selection gradient level
ssGradLevel = grSliceIR(round(mean(rfLimitsIR)));

%% verify feasibility
% minimum TI is from center of RF to end to RF 
timeCenterRF = 0.5*rfIRTime(rfLimitsIR(2));
if rfLimitsIR > 1
     timeCenterRF = timeCenterRF + rfIRTime(rfLimitsIR(1)-1)*0.5;
end
minTI = rfIRTime(end) - timeCenterRF;
if prepIR.TI <= 0
    % use minimum possible TI
    prepIR.TI = minTI;
end
    tWaitIR  = prepIR.TI - minTI;
if tWaitIR < 0
    msg = sprintf( ['The selected Inversion Recovery Time (TI=%.2fms) is too short. ',...
        'Minimum Inversion Recovery Time for the current configuration is TI=%.2fms'],...
        prepIR.TI*1e3, minTI*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
else
    tWaitIR = max(tWaitIR,1e-12);
end


%% combine the into a sequence
numSteps    = length(rfIRTime);

% start/end of RF IR
iniRFIR         = 1;
endRFIR         = numSteps;
% total time of repetition
numSteps        = numSteps + 1;

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
repetition.time(iniRFIR:endRFIR)     = rfIRTime;
%% RF
repetition.rfEntries(rfLimitsIR(1):rfLimitsIR(2))   = 1;
repetition.rfmSignal(iniRFIR:endRFIR)  	= rfmIR;
repetition.rfpSignal(iniRFIR:endRFIR)  	= rfpIR;
repetition.rffSignal(iniRFIR:endRFIR)  	= rffIR;
repetition.ssSignal(iniRFIR:endRFIR)   	= grSliceIR;

%% assign end time
repetition.time(end)                    = rfIRTime(end) + tWaitIR;


%% pulse sequences info
repetition.type         = 'prep-IR';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = ssGradLevel;

%% for splitting into parts (RF vs Non-RF)
% multiple parts, depending on its characteristics
numParts = 1;
if rfLimitsIR(1) > 1
    numParts = numParts + 1;
    partLimits = [1, rfLimitsIR(1)-1; rfLimitsIR];
    partType{1} = 'GR';
    partType{2} = 'RF';
else
    partLimits = rfLimitsIR;
    partType{1} = 'RF';
end
% first gradient
numParts = numParts + 1;
partLimits = [partLimits; rfLimitsIR(2)+1, repetition.numSteps];
%% apply final SWC if required
if prepIR.postRFswc
    % if we have SWC we need RO part to be taken into account
    repetition.swcSignal(end)   = 1;
    partType{numParts}          = 'RO';
else
    partType{numParts}          = 'GR';
end

% assign
repetition.numParts   = numParts; % number of parts
repetition.partType   = partType; % type of part
repetition.partLimits = partLimits; % index start/end of part
repetition.rxLimits   = []; % index start/end of readouts
%
repetition.TI = prepIR.TI;

%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s repetition, elapsed time %.3fs',...
        functionName, repetition.type, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', repetition.totalTime);
    fprintf(fid, '\n  Number steps        %d', repetition.numSteps);
    fprintf(fid, '\n  Number readouts     %d', repetition.numRxs);
    fprintf(fid, '\n  Effective TI        %.3fms', 1e3*repetition.TI );
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


