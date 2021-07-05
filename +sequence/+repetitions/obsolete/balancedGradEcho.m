function [repetition] = balancedGradEcho(...
    acqData, mainRF, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.BALANCEDDGRADECHO
%
%	Generates a Gradient Echo repetition, with fully-balanced encoding.
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
functionName = 'sequence.repetitions.balancedGradEcho';

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
[encTime,signalFE,signalPE,signalSE,signalRX,dtEnc,effTE,rxLimits] = ...
    sequence.encodings.cartesianFullyBalanced(...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
    acqData.numFE, acqData.numPE, acqData.numSE,...
    acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
    expControl.sequence.dtGR, expControl );

%% generate the waveforms for the RF (use main RF info)
[rfTime,rfMag,rfPhase,grSlice,rfLimits] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% get slice selection gradient level
ssGradLevel = grSlice(round(mean(rfLimits)));

%% verify feasibility of TE/TR and correctness of timing
% TE
minTE = rfTime(end) - mean(rfTime(rfLimits)) + effTE;
if acqData.forceMinTE
    % force to minimum TE possible
    tWaitRF = 0;
    acqData.TE = minTE;
else
    tWaitRF  = acqData.TE - minTE;
end
% TR
minTR = rfTime(end) + encTime(end) + tWaitRF;
if acqData.forceMinTR
    % force to minimum TR possible
    tWaitRep = 0;
    acqData.TR = minTR;
else
    tWaitRep = acqData.TR - minTR;
end
if tWaitRF < 0
    if ~isempty(expControl.connLocalDB) ...
            && strcmpi(expControl.application,'edutool')
        msg = sprintf( ['The selected Echo Time (TE=%.3fms) is too short. ',...
            'Minimum Echo Time for the current configuration is TE=%.3fms'],...
            acqData.TE*1e3, minTE*1e3 );
        eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error', expControl.experimentID, expControl.pulseqID);
    else
        ME = MException('sequence:wrongEchoTime',...
            '%s : %.3fms TE too short (minimum %.3fms)',...
            functionName, acqData.TE*1e3, minTE*1e3);
        throw(ME);
    end
else
    tWaitRF = max(tWaitRF,1e-12);
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
numTimesRF    = length(rfTime);
numTimesEnc   = length(encTime);

% start/end of RF90
iniRF           = 1;
endRF           = iniRF + numTimesRF-1;
rfLimits        = rfLimits + iniRF-1;
% start/end of Encoding gradients
idxWaitRF       = endRF + 1;
iniEnc          = idxWaitRF + 1;
endEnc          = iniEnc+numTimesEnc-1;
idxWaitRep      = endEnc + 1;
% correct rxLimits
rxLimits        = rxLimits + idxWaitRF;
% total time of repetition
numSteps        = idxWaitRep;

%% repetition waveforms
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


%% assign blocks
repetition.time(iniRF:endRF)         = rfTime;
repetition.time(idxWaitRF)           = tWaitRF + repetition.time(idxWaitRF-1);
repetition.time(iniEnc:endEnc)       = encTime + repetition.time(iniEnc-1);
repetition.time(idxWaitRep)          = tWaitRep + repetition.time(idxWaitRep-1);
%% RF
repetition.rfEntries(rfLimits(1):rfLimits(2))   = 1;
repetition.rfmSignal(iniRF:endRF)    = rfMag;
repetition.rfpSignal(iniRF:endRF)    = rfPhase;
repetition.ssSignal(iniRF:endRF)     = grSlice;
%% Encoding
repetition.feSignal(iniEnc:endEnc) = ...
    repetition.feSignal(iniEnc:endEnc) + signalFE(:);
repetition.peSignal(iniEnc:endEnc) = ...
    repetition.peSignal(iniEnc:endEnc) + signalPE(:);
repetition.seSignal(iniEnc:endEnc) = ...
    repetition.seSignal(iniEnc:endEnc) + signalSE(:);
repetition.rxSignal(iniEnc:endEnc) = signalRX;

%% pulse sequences info
repetition.type         = 'Balanced Gradient Echo';
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
    partLimits = [1, rfLimits(2)];
    partType{numParts} = 'RF';
else
    % multiple parts, depending on its characteristics
    numParts = 1;
    % RF
    if rfLimits(1) > 1
        numParts = numParts + 1;
        partLimits = [1, rfLimits(1)-1; rfLimits];
        partType{1} = 'GR';
        partType{2} = 'RF';
    else
        partLimits = rfLimits;
        partType{1} = 'RF';
    end
end
% rest may have readouts
numParts = numParts + 1;
partLimits = [partLimits; rfLimits(2)+1, repetition.numSteps];
partType{numParts} = 'RO';
% assign
repetition.numParts   = numParts; % number of parts
repetition.partType   = partType; % type of part
repetition.partLimits = partLimits; % index start/end of part
repetition.rxLimits   = rxLimits; % index start/end of readouts
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
%         plot(repetition.time(rfLimits),[0,0], '^');
%         plot(repetition.time,repetition.feSignal);
%         plot(repetition.time,repetition.peSignal);
%         plot(repetition.time,repetition.ssSignal);
%         plot(repetition.time(repetition.rxSignal>0), ...
%             repetition.feSignal(repetition.rxSignal>0), 'o')

end


