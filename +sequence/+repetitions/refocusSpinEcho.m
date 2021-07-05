function [repetition] = refocusSpinEcho(...
    acqData, refRF, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.REFOCUSSPINECHO
%
%	Generates refocusing Encoding for Spin Echo family.
%   RF180 +PE  +FE  -PE 
%
% INPUT
%   acqData
%   refRF
%   mrSystem        
%   expControl      
%
% OUTPUT
%   repetition   repetition struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:repetitions:refocusSpinEcho';

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
% [encTime,signalFE,signalPE,signalSE,signalRX,areaFE,dtEnc,effTE,rxLimits] = ...
%     sequence.encodings.cartesianPhaseBalanced( ...
%     acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
%     acqData.numFE, acqData.numPE, acqData.numSE,...
%     acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
%     mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
%     expControl.sequence.dtGR, expControl );

% define vars
freqBalance         = 'none';
phaseBalance        = 'full';
sliceBalance        = 'full';
echoTimes           = []; 
reverseEchoPolarity = 1;
% generate encoding using general function
[encTime,signalFE,signalPE,signalSE,signalRX,rxLimits,effTE,areaFE,dtEnc] = ...
    sequence.encodings.cartesian( ...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
    acqData.numFE, acqData.numPE, acqData.numSE,...
    acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
    freqBalance, phaseBalance, sliceBalance, ...
    echoTimes, reverseEchoPolarity,... 
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
    expControl.sequence.dtGR, expControl );


%% generate the waveforms for the RF 180 (use refocusing RF info)
[rfTime,rfm,rfp,rff,grSlice,rfLimits] = ...
    sequence.excitations.sliceSelectRF( refRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% get slice selection gradient level in the 180 
ssGradLevel = grSlice(round(mean(rfLimits)));

%% verify feasibility of TE/TR and correctness of timing
% min TE is the maximum on
%   time between pulses + preRW : (estimated) rfTime(end) + timeHalfFE
%   time between center of 180 and effective time

% RF guard as the maximum of user defined and the minimum required
expControl.sequence.tGuardRF = 10e-6;
tGuardRF = max(expControl.sequence.tGuardRF, 1e-12);
tEndRF   = max(rfTime(rfLimits(:))); % end of RF pulse

% % Estimate minimum TE for the RF90
% %    time to generate half the FE area at half stregth (preRW)
% timeHalfFE = areaFE/mrSystem.maxGStrenght;
% minTE90    = tEndRF + tGuardRF + timeHalfFE;

% get minimum TE for the 180
tStartRO = min(encTime(rxLimits(:))) - dtEnc; 
% find min guard time so that SS gradient do not overlap RO samples
tSSRW    = rfTime(end) - tEndRF; % duration of Slice selection rewinder after RF
tGuardRF = max((tSSRW - tStartRO), tGuardRF);
% update tEndRF to min time at which encoding can start (add guard)
tEndRF   =  tEndRF + tGuardRF;
% minimun TE possible (so that there is no overlap with RF)
minTE180 = tEndRF - mean(rfTime(rfLimits)) + effTE;

% minimum TE 
minTE = minTE180;
if acqData.forceMinTE || acqData.TE/2 < minTE
    % force to minimum TE possible
    tWaitRF = 0;
    repTE   = 2*minTE; % actual repetition TE
else
    tWaitRF = acqData.TE/2 - minTE;
    repTE   = acqData.TE; % actual repetition TE
end

% minimum TR possible, including wait to get TE
minTR = tEndRF + tWaitRF + encTime(end);
if repTE < minTR 
    % TR is less than minimum
    tWaitRep    = 0.0;
    tWaitRF     = tWaitRF + (minTR - repTE)/2; % update echo wait
    repTE = minTR; % update TE to actual
else
    tWaitRep = repTE - minTR;
end

% add small value if 0 wait
tWaitRF = max(tWaitRF,1e-12);
tWaitRep = max(tWaitRep,1e-12);

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
    % time points after the end of the encoding:
    iQuery = commonIndex >= encIndex(end);
    %   zero them
    repetition.feSignal(commonIndex(iQuery)) = 0.0;
    repetition.peSignal(commonIndex(iQuery)) = 0.0;
    repetition.seSignal(commonIndex(iQuery)) = 0.0;
    
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
repetition.rfmSignal(rfIndex)   = rfm;
repetition.rfpSignal(rfIndex)   = rfp;
repetition.rffSignal(rfIndex)   = rff;
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
repetition.type         = 'refocus-Spin-Echo';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = ssGradLevel;

%% for splitting into parts (RF vs Non-RF)
% multiple parts, depending on its characteristics
numParts = 1;
% 180
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
repetition.TE         = repTE;
repetition.TR         = repetition.time(end);

%% verify PE area
if abs(sum(repetition.timeDiff.*repetition.peSignal)) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : PE gradient area mismatch in Spin Echo encoding',functionName);
    throw(ME);
end
%% verify FE area
if abs( sum(repetition.timeDiff.*repetition.feSignal) - areaFE ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : FE gradient area mismatch in Spin Echo encoding',functionName);
    throw(ME);
end
%% verify SE area
if abs( sum(repetition.timeDiff.*repetition.seSignal) ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : SE gradient area mismatch in Spin Echo encoding',functionName);
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
T180toEcho = mean(repetition.time(repetition.rxSignal>0)) ...
    - mean(repetition.time(rfLimits));
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


