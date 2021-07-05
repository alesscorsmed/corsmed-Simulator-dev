function [repetition] = unbalancedFullSpinEcho(...
    acqData, mainRF, refRF, mrSystem, expControl)
%
% SEQUENCE.REPETITIONS.UNBALANCEDFULLSPINECHO
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
functionName = 'sequence:repetitions:unbalancedFullSpinEcho';

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

%% generate general encodings
% define vars
freqBalance         = 'none';
phaseBalance        = 'none';
sliceBalance        = 'none';
echoTimes           = [];
reverseEchoPolarity = 1;
% generate encoding using general function
[roTime, roSignalFE, ~,~, roSignalRX, rxLimits, effTE, areaFE, dtEnc] = ...
    sequence.encodings.cartesian( ...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, acqData.fovSE,...
    acqData.numFE, acqData.numPE, acqData.numSE,...
    acqData.samplingFactorFE, acqData.samplingFactorPE, acqData.samplingFactorSE,...
    freqBalance, phaseBalance, sliceBalance, ...
    echoTimes, reverseEchoPolarity,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma,...
    expControl.sequence.dtGR, expControl );

%% generate the rewinders
areaPE = acqData.numPE/(acqData.gamma*acqData.fovPE)/2;
areaSE = acqData.numSE/(acqData.gamma*acqData.fovSE)/2;
[ rwTime,rwSignalFE,rwSignalPE,rwSignalSE ] = ...
    sequence.encodings.cartesianRewinder( areaFE/2, areaPE, areaSE, ...
    mrSystem.maxGStrenght, mrSystem.SlewRate, dtEnc, expControl);

%% generate the waveforms for the RF 90 (use main RF info)
[rf90Time,rfm90,rfp90,rff90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% generate the waveforms for the RF 180 (use refocusing RF info)
refRF.sliceThickness = refRF.sliceThickness*3;
[rf180Time,rfm180,rfp180,rff180,grSlice180,rfLimits180] = ...
    sequence.excitations.sliceSelectRF( refRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% get slice selection gradient level in the 180 
ssGradLevel = grSlice180(round(mean(rfLimits180)));

%% RF guard as the maximum of user defined and the minimum required
expControl.sequence.tGuardRF = 1e-6;
tGuardRF = max(expControl.sequence.tGuardRF, 1e-12);

%% verify feasibility of TE and correctness of timing for Refocusing
% get minimum TE for the 180
tEnd180  = max(rf180Time(rfLimits180(:))); % end of RF pulse
tStartRO = min(roTime(rxLimits(:))) - dtEnc; 
% find min guard time so that SS gradient do not overlap RO samples
tSSRW    = rf180Time(end) - tEnd180; % duration of Slice selection rewinder after RF
tGuardRF = max((tSSRW - tStartRO), tGuardRF);
% update tEndRF to min time at which encoding can start (add guard)
tEnd180  =  tEnd180 + tGuardRF;
% minimun TE possible (so that there is no overlap with RF)
minTE180 = tEnd180 - mean(rf180Time(rfLimits180)) + effTE;

%% verify feasibility of TE and correctness of timing for 90
tEnd90   = max(rf90Time(rfLimits90(:))); % end of RF pulse
% update tEndRF to min time at which encoding can start (add guard)
tGuardRF = max(expControl.sequence.tGuardRF, 1e-12);
tEnd90   =  tEnd90 + tGuardRF;
% minimun TE possible (so that there is no overlap with RF)
minTR90  = max(tEnd90 + rwTime(end), rf90Time(end));
minTE90  = minTR90 - mean(rf90Time(rfLimits90)) + mean(rf180Time(rfLimits180));

minTE = 2*max(minTE180,minTE90);
if acqData.forceMinTE
    % force to minimum TE possible
    acqData.TE = minTE;
end
tWait90  = acqData.TE/2 - minTE90;
tWait180 = acqData.TE/2 - minTE180;

if acqData.TE < minTE
    msg = sprintf( ['The selected Echo Time (TE=%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is TE=%.2fms'],...
        acqData.TE*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end
    
tWait180 = max(tWait180,1e-12);
tWait90 = max(tWait90,1e-12);

%% generate the pre encoding repetition
shift90 = tEnd90 + tWait90/2;
repTime90 = minTR90 + tWait90;
[repetition90] = sequence.tools.combineIntoRepetition( ...
    rf90Time, rfm90, rfp90, rff90, grSlice90, rfLimits90, ...
    rwTime, rwSignalFE, -rwSignalPE, -rwSignalSE, 0*rwSignalSE, [], ...
    shift90, repTime90 );

%% generate the encoding repetition
shift180   = tEnd180 + tWait180;
repTime180 = acqData.TE + 1e-12;
[repetition180] = sequence.tools.combineIntoRepetition( ...
    rf180Time, rfm180, rfp180, rff180, grSlice180, rfLimits180, ...
    roTime, roSignalFE, 0*roSignalFE, 0*roSignalFE, roSignalRX, rxLimits, ...
    shift180, repTime180 );

%% verify feasibility of TR and correctness of timing
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

%% verify SS area
if abs( sum(repetition.timeDiff.*repetition.ssSignal) - ...
        sum(reshape( diff([0.0; rf90Time]),[],1 ).*grSlice90) - ...
        sum(reshape( diff([0.0; rf180Time]),[],1 ).*grSlice180) ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : Slice selection gradient area mismatch',functionName);
    throw(ME);
end
%% verify GR areas
if abs( sum(repetition.timeDiff.*repetition.seSignal) + ...
        sum(reshape( diff([0.0; rwTime]),[],1 ).*rwSignalSE) ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : Slice encoding gradient area mismatch',functionName);
    throw(ME);
end
if abs( sum(repetition.timeDiff.*repetition.peSignal) + ...
        sum(reshape( diff([0.0; rwTime]),[],1 ).*rwSignalPE) ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : Phase encoding gradient area mismatch',functionName);
    throw(ME);
end
if abs( sum(repetition.timeDiff.*repetition.feSignal) - 3*areaFE/2 ) > 1e-12
    ME = MException('sequence:wrongGradArea',...
        '%s : Frequency encoding gradient area mismatch',functionName);
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
%         plot(repetition.time(repetition.partLimits(:,1)),0*repetition.partLimits(:,1), '^');
%         plot(repetition.time(repetition.partLimits(:,2)),0*repetition.partLimits(:,2), 'v');
%         plot(repetition.time,repetition.feSignal);
%         plot(repetition.time,repetition.peSignal);
%         plot(repetition.time,repetition.ssSignal);
%         plot(repetition.time(repetition.rxSignal>0), ...
%             repetition.feSignal(repetition.rxSignal>0), 'o')

end


