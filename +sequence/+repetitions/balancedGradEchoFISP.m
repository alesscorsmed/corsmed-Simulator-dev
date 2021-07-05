function [repetition] = balancedGradEchoFISP(...
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
functionName = 'sequence:repetitions:balancedGradEcho';

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
% [encTime,signalFE,signalPE,signalSE,signalRX,dtEnc,effTE,rxLimits] = ...
%     sequence.encodings.cartesianFullyBalanced(...
%     acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
%     acqData.numFE, acqData.numPE, acqData.numSE,...
%     acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
%     mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
%     expControl.sequence.dtGR, expControl );

% define vars
freqBalance         = 'spoil';
phaseBalance        = 'full';
sliceBalance        = 'full';
echoTimes           = []; 
reverseEchoPolarity = 1;
% generate encoding using general function
[encTime,signalFE,signalPE,signalSE,signalRX,rxLimits,effTE,~,dtEnc] = ...
    sequence.encodings.cartesian( ...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
    acqData.numFE, acqData.numPE, acqData.numSE,...
    acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
    freqBalance, phaseBalance, sliceBalance, ...
    echoTimes, reverseEchoPolarity,... 
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
    expControl.sequence.dtGR, expControl );


%% generate the waveforms for the RF (use main RF info)
[rfTime,rfMag,rfPhase,rfFreq,grSlice,rfLimits] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% get slice selection gradient level
ssGradLevel = grSlice(round(mean(rfLimits)));

%% verify feasibility of TE/TR and correctness of timing
% time of sampling starts in Encoding part
tStartRO = min(encTime(rxLimits(:))) - dtEnc; 
% find min guard time so that SS gradient do not overlap RO samples
tEndRF   = max(rfTime(rfLimits(:))); % end of RF pulse
tSSRW    = rfTime(end) - tEndRF; % duration of Slice selection rewinder after RF
tGuardRF = max((tSSRW - tStartRO), 1e-12);
% RF guard as the maximum of user defined and the minimum required
expControl.sequence.tGuardRF = 10e-6;
tGuardRF = max(expControl.sequence.tGuardRF, tGuardRF);
% update tEndRF to min time at which encoding can start (add guard)
tEndRF   =  tEndRF + tGuardRF;
% minimun TE possible (so that there is no overlap with RF)
minTE    = tEndRF - mean(rfTime(rfLimits)) + effTE;
if acqData.forceMinTE
    % force to minimum TE possible
    tWaitRF = tGuardRF;
    acqData.TE = minTE;
else
    % calculate the wait to match the required TE
    tWaitRF  = acqData.TE - minTE;
end
% minimum TR possible, including wait to get TE
minTR = tEndRF + tWaitRF + encTime(end);
if acqData.forceMinTR
    % force to minimum TR possible
    tWaitRep = 0;
    acqData.TR = minTR;
else
    tWaitRep = acqData.TR - minTR;
end
if tWaitRF < 0
    msg = sprintf( ['The selected Echo Time (TE=%.3fms) is too short. ',...
        'Minimum Echo Time for the current configuration is TE=%.3fms'],...
        acqData.TE*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
else
    tWaitRF = max(tWaitRF,1e-12);
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
else
    tWaitRep = max(tWaitRep,1e-12);
end

%% correct encoding times to generate correct TE
midTime = tEndRF + tWaitRF;
encTime = encTime + midTime;

%% append midtime as part of the encoding, with zero signals
encTime   = [midTime; encTime(:)];
signalFE = [0.0; signalFE(:)];
signalPE = [0.0; signalPE(:)];
signalSE = [0.0; signalSE(:)];
signalRX = [0.0; signalRX(:)];
rxLimits = rxLimits + 1; % added an extra point at beginning

%% combine into the repetition: unify times
repetition.time = reshape(union([rfTime; midTime], encTime),[],1);
% increase the number to match the required TR
numSteps = length(repetition.time) + 1;
repetition.time(numSteps,1) = repetition.time(end) + tWaitRep;
% time differences
repetition.timeDiff = reshape(diff([0; repetition.time]),[],1);

%% allocate waveforms
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

%% find incidence of RF and ENC times into new time array
% get indexes of RF arrays into new time array
[~,~,rfIndex] = intersect(rfTime,repetition.time,'stable');
% get indexes of RF arrays into new time array
[~,~,encIndex] = intersect(encTime,repetition.time,'stable');

%% find intertwinned entries (if any)
commonIndex = encIndex(1):rfIndex(end);
% assign overlapping signals: make sure to keep area
if ~isempty(commonIndex)
    
    %% interpolate the Encoding signals
    % time points before the starting of the encoding:
    iQuery = commonIndex <= encIndex(1);
    %   zero them
    repetition.feSignal(commonIndex(iQuery)) = 0.0;
    repetition.peSignal(commonIndex(iQuery)) = 0.0;
    repetition.seSignal(commonIndex(iQuery)) = 0.0;
    % time points after the starting of the encoding
    iQuery = commonIndex > encIndex(1);
    %   interpolate using next neighbor: keep areas 
    repetition.feSignal(commonIndex(iQuery)) = interp1( encTime, signalFE,...
        repetition.time(commonIndex(iQuery)), 'next');
    repetition.peSignal(commonIndex(iQuery)) = interp1( encTime, signalPE,...
        repetition.time(commonIndex(iQuery)), 'next');
    repetition.seSignal(commonIndex(iQuery)) = interp1( encTime, signalSE,...
        repetition.time(commonIndex(iQuery)), 'next');
    
    %% interpolate the SS gradient points
    % time points before the end of the ss:
    iQuery = commonIndex <= rfIndex(end);
    %   interpolate using next neighbor: keep areas 
    repetition.ssSignal(commonIndex(iQuery)) = interp1( rfTime, grSlice,...
        repetition.time(commonIndex(iQuery)), 'next');
    % time points after the end of the ss:
    iQuery = commonIndex > rfIndex(end);
    %   zero them
    repetition.ssSignal(commonIndex(iQuery)) = 0.0;

end

%% assign RF signals
% assign signal entries to corrsponging positions
repetition.rfmSignal(rfIndex)   = rfMag;
repetition.rfpSignal(rfIndex)   = rfPhase;
repetition.rffSignal(rfIndex)   = rfFreq;
repetition.ssSignal(rfIndex)    = grSlice;
% correct the rfLimits (if needed) and define the actual RF entries
rfLimits(:) = rfIndex(rfLimits(:));
repetition.rfEntries(rfLimits(1):rfLimits(2)) = 1;

%% assign ENC signals
% assign signal entries to corrsponging positions
repetition.feSignal(encIndex) = signalFE(:);
repetition.peSignal(encIndex) = signalPE(:);
repetition.seSignal(encIndex) = signalSE(:);
repetition.rxSignal(encIndex) = signalRX(:);
% correct the position of the rxLimits
rxLimits(:) = encIndex(rxLimits(:));

%% pulse sequences info
repetition.type         = 'FISP';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = ssGradLevel;

%% for splitting into parts (RF vs Non-RF)
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
% add info of minimum TR (including waiting for correct TE)
repetition.minTR      = minTR;

%% verify PE area
if abs(sum(repetition.timeDiff.*repetition.peSignal)) > 1e-12 && lower(phaseBalance) ~= "spoil"
    ME = MException('sequence:wrongGradArea',...
        '%s : PE gradient area mismatch in Balanced GRE encoding',functionName);
    throw(ME);
end
%% verify FE area
if abs(sum(repetition.timeDiff.*repetition.feSignal)) > 1e-12 && lower(freqBalance) ~= "spoil"
    ME = MException('sequence:wrongGradArea',...
        '%s : FE gradient area mismatch in Balanced GRE encoding',functionName);
    throw(ME);
end
%% verify SE area
if abs(sum(repetition.timeDiff.*repetition.seSignal)) > 1e-12 && lower(sliceBalance) ~= "spoil"
    ME = MException('sequence:wrongGradArea',...
        '%s : SE gradient area mismatch in Balanced GRE encoding',functionName);
    throw(ME);
end
%% verify SS area
if abs( sum(repetition.timeDiff.*repetition.ssSignal) - ...
        sum(reshape( diff([0.0; rfTime]),[],1 ).*grSlice) ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : Slice selection gradient area mismatch',functionName);
    throw(ME);
end
%% verify echo times
TRFtoEcho = mean(repetition.time(repetition.rxLimits(1,:))) ...
    - mean(repetition.time(rfLimits));
if abs(TRFtoEcho - repetition.TE ) > 1e-12
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
%         figure(1);
%         plot(repetition.time,repetition.rfmSignal*1e3)
%         hold on
%         plot(repetition.time(repetition.rfEntries),...
%             0*repetition.rfEntries, '^');
%         plot(repetition.time,repetition.feSignal);
%         plot(repetition.time,repetition.peSignal);
%         plot(repetition.time,repetition.ssSignal);
%         plot(repetition.time(repetition.rxSignal>0), ...
%             repetition.feSignal(repetition.rxSignal>0), 'o')


end


