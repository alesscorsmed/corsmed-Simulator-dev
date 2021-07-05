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
functionName = 'sequence:repetitions:prepSpinEcho';

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
[rfTime,rfm,rfp,rff,grSlice,rfLimits] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl);

%% get slice selection gradient level in the 180 
ssGradLevel = grSlice(round(mean(rfLimits)));

%% verify feasibility of TE and correctness of timing
% RF guard as the maximum of user defined and the minimum required
expControl.sequence.tGuardRF = 10e-6;
tGuardRF = max(expControl.sequence.tGuardRF, 1e-12);
% find min guard time so that SS gradient do not overlap RO samples
tEndRF   = max(rfTime(rfLimits(:))); % end of RF pulse
% update tEndRF to min time at which encoding can start (add guard)
tEndRF   =  tEndRF + tGuardRF;
% minimun TE possible (so that there is no overlap with RF)
minTE90  = max(tEndRF + rwTime(end), rfTime(end));

if acqData.TE/2 < minTE90
    % force to minimum TE possible
    tWait90 = 0;
    repTE   = 2*minTE90; % actual repetition TE
else
    tWait90  = acqData.TE/2 - minTE90;
    repTE   = acqData.TE; % actual repetition TE
end
tWait90 = max(tWait90,1e-12);

%% divide the twait into after RF and after the RW
tWaitRF  = tWait90/2;
tWaitRep = tWait90/2;

%% correct RW times to generate correct TE
midTime = tEndRF + tWaitRF;
rwTime = rwTime + midTime;

%% append midtime as part of the encoding, with zero signals
rwTime   = [midTime; rwTime(:)];
rwSignal = [0.0; rwSignal(:)];

%% combine into the repetition: unify times
repetition.time = reshape(union(rfTime, rwTime),[],1);
% increase the number to match the required TR
numSteps = length(repetition.time) + 1;
repetition.time(numSteps,1) = repTE/2;
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
[~,~,rwIndex] = intersect(rwTime,repetition.time,'stable');

%% find intertwinned entries (if any)
commonIndex = rwIndex(1):rfIndex(end);
% assign overlapping signals: make sure to keep area
if ~isempty(commonIndex)
    
    %% interpolate the Encoding signals
    % time points before the starting of the encoding:
    iQuery = commonIndex <= rwIndex(1);
    %   zero them
    repetition.feSignal(commonIndex(iQuery)) = 0.0;
    % time points after the starting of the encoding
    iQuery = commonIndex > rwIndex(1);
    %   interpolate using next neighbor: keep areas 
    repetition.feSignal(commonIndex(iQuery)) = interp1( rwTime, rwSignal,...
        repetition.time(commonIndex(iQuery)), 'next');
    % time points after the end of the encoding:
    iQuery = commonIndex >= rwIndex(end);
    %   zero them
    repetition.feSignal(commonIndex(iQuery)) = 0.0;
    
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

%% assign RW signals
% assign signal entries to corrsponging positions
repetition.feSignal(rwIndex) = rwSignal(:);

%% pulse sequences info
repetition.type         = 'prep-Spin-Echo';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = ssGradLevel;

%% for splitting into parts (RF vs Non-RF)
% multiple parts, depending on its characteristics
numParts = 1;
% 90
if rfLimits(1) > 1
    numParts = numParts + 1;
    partLimits = [1, rfLimits(1)-1; rfLimits];
    partType{1} = 'GR';
    partType{2} = 'RF';
else
    partLimits = rfLimits;
    partType{1} = 'RF';
end
% first gradient
numParts = numParts + 1;
partLimits = [partLimits; rfLimits(2)+1, repetition.numSteps];
partType{numParts} = 'GR';
% assign
repetition.numParts   = numParts; % number of parts
repetition.partType   = partType; % type of part
repetition.partLimits = partLimits; % index start/end of part
repetition.rxLimits   = []; % index start/end of readouts
% update TE / TR
repetition.TE         = repTE;
repetition.TR         = repetition.time(end);

%% verify FE area
sum(reshape( diff([0.0; rwTime]),[],1 ).*rwSignal);
if abs(sum(repetition.timeDiff.*repetition.feSignal)-areaFE/2) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : FE gradient area mismatch in Spin Echo FE Rewinder',functionName);
    throw(ME);
end
%% verify SS area
if abs( sum(repetition.timeDiff.*repetition.ssSignal) - ...
        sum(reshape( diff([0.0; rfTime]),[],1 ).*grSlice) ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : Slice selection gradient area mismatch',functionName);
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


