function [pulseSequence] = MOLLI(...
    acquisition, encoding, mrSystem, expControl, anatomicalModel)
%
% SEQUENCE.FAMILYSSFP.MOLLI
%
%	Generates a MOLLI sequence.
%
% INPUT
%   acquisition         
%   mrSystem        
%   expControl      
%
% OUTPUT
%   pulseSequence   pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.familySSFP.MOLLI';

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

%% get required info
acqData     = acquisition.data;
prepIR      = acquisition.prepIR;

%% extract info for MOLLI
cardiacCycleDurationMsec = (60/anatomicalModel.HR)*1e3;

% compute the time for each cardiac cycle (duration of the bSSFP shot)
cardiacShotTR = round(cardiacCycleDurationMsec*1e-3); % cardiac cycle duration from ms to s

% molli scheme
[schemeMOLLI,pauseMOLLI,TIs] = ...
    sequence.tools.extractMOLLIspecs(acquisition.data.MOLLI,...
    cardiacCycleDurationMsec);

% modify a little the molliScheme
molliShots  = 0;
molliNumTI  = 0;
for rr = 1:numel(schemeMOLLI)
    molliScheme{rr}.shots   = schemeMOLLI(rr);
    molliScheme{rr}.pause   = pauseMOLLI;
    molliScheme{rr}.TI      = TIs(rr);
    molliShots              = molliShots + schemeMOLLI(rr);
    if TIs(rr) > 0
        molliNumTI          = molliNumTI + 1;
    end
end

%% create the bSSFP building block
[bssfpSequence] = sequence.familySSFP.balancedSSFP(...
    acquisition, encoding, mrSystem, expControl);

% check if the number of reps is even to modify polarity
evenBssfpReps = rem(bssfpSequence.numReps+1,2);

%% create the IR with minimum TI
prepIR.TI       = -1; % this will generate the minimum TR
prepIR.postRFswc=  1; % apply SWC after RF
repetitionIR    = sequence.repetitions.modPrepInvRecovery(...
    acqData, prepIR, mrSystem, expControl);

%% compute the minimum TI: from end of RF to center of bSSFP
tEndIR  = max(repetitionIR.time(repetitionIR.rfEntries > 0)); 
minTI   = tEndIR + bssfpSequence.effTE;

%% allocate space for the full sequence
[pulseSequence] = data.simulation.initializeSequence();

%% NOTE: we add one extra point per bssfp sequence
%% this extra point will have zero entries, and its timing will be 
%% just a small delta before the start of the next bssfp (or IR)
%% This avoids cases in which we may have large RF steps 
%% (if the first point of the RF is not pure gradient)
% total number of steps
numSteps = molliShots*bssfpSequence.numSteps ...
    + molliNumTI*repetitionIR.numSteps + molliShots;
% total number of parts
numParts = molliShots*bssfpSequence.numParts ...
    + molliNumTI*repetitionIR.numParts;
% steps
pulseSequence.numSteps  = numSteps;
% initialize data
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.gxSignal  = zeros(numSteps,1);
pulseSequence.gySignal  = zeros(numSteps,1);
pulseSequence.gzSignal  = zeros(numSteps,1);
pulseSequence.rfmSignal = zeros(numSteps,1);
pulseSequence.rfpSignal = zeros(numSteps,1);
pulseSequence.rffSignal = zeros(numSteps,1);
pulseSequence.swcSignal = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1);
% parts
pulseSequence.numParts = numParts;
pulseSequence.partType{pulseSequence.numParts} = []; % allocate
pulseSequence.partLimits = zeros(pulseSequence.numParts,2);
% rx Limits: use number of encodings w/ readouts
pulseSequence.rxLimits = zeros(molliShots*bssfpSequence.numEnc,2);

%% loop on the molli scheme, and assemble the pulseSequence
currentShot     = 0;
currentCycle    = 0;
currentTime     = 0;
triggerTime     = 0;
endShotTime     = 0;
% initiate
actualTI        = zeros(1,sum(schemeMOLLI));
% cumulative shifts
idxShift  = 0;
partShift = 0;
limitShift= 0;
% to it
for rr = 1:numel(molliScheme)
    % shot time
    shotTime = bssfpSequence.totalTime;
    %% verify timings
    if molliScheme{rr}.TI > 0
        % in case of IR
        tWaitIR = molliScheme{rr}.TI - minTI;
        % verify correct timing
        if tWaitIR < 0
            msg = sprintf( ['The selected IR time (TI=%.2fms) is too short. ',...
                'Minimum TI for the current configuration is TI=%.2fms'],...
                molliScheme{rr}.TI*1e3, minTI*1e3 );
            if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
                eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
            else
                ME = MException('error:user', '%s', msg);
                throw(ME);
            end
        end
        % correct shot time with IR part
        shotTime = shotTime + repetitionIR.totalTime + tWaitIR;
    end
    % verify total time of shot
    if shotTime > cardiacShotTR
        msg = sprintf( ['The shot duration (%.2fms) is too short. ',...
            'Cardiac cycle duration is %.2fms'],...
            shotTime*1e3, cardiacShotTR*1e3 );
        if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
            eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
        else
            ME = MException('error:user', '%s', msg);
            throw(ME);
        end
    end
    
    %% all peachy: add shot parts
    
    % add the IR if needed
    if molliScheme{rr}.TI > 0
        % target start time for shot
        targetStartTime = triggerTime + currentCycle*cardiacShotTR;
        % get total time of IR block:
        irTime = repetitionIR.time(end) + tWaitIR;
        % if it is not the first IR, correct current time
        % so that bSSFP starts at multiple of cardiacShotTR
        if rr > 1
            targetStartTime = targetStartTime - irTime;
            % set as time for the current index (last point of prev bssfp)
            % a time just before starting the IR: 
            % avoid potential large RF steps
            pulseSequence.time(idxShift) = targetStartTime - 1e-8; 
        end
        % check correct timings
        if endShotTime > targetStartTime
            % Throw error: we are computing things badly
            msg = sprintf( ['(in IR block) Time end of prev. bssfp (%.2fms)',...
                ' larger than target time (%.2fms) for cycle %d'],...
                endShotTime*1e3, targetStartTime*1e3, currentCycle );
            ME = MException('error:edutool', '%s', msg);
            throw(ME);
        end
         % starting time at the desired time
        currentTime = targetStartTime;
        % target indexes
        tgtIdx  = (1:repetitionIR.numSteps) + idxShift;
        tgtPart = (1:repetitionIR.numParts) + partShift;
        % assign waveforms
        % gradients
        pulseSequence.gxSignal(tgtIdx)      = repetitionIR.feSignal;
        pulseSequence.gySignal(tgtIdx)      = repetitionIR.peSignal;
        pulseSequence.gzSignal(tgtIdx)      = repetitionIR.seSignal;
        % rf
        pulseSequence.rfmSignal(tgtIdx)     = repetitionIR.rfmSignal;
        pulseSequence.rfpSignal(tgtIdx)     = repetitionIR.rfpSignal;
        pulseSequence.rffSignal(tgtIdx)     = repetitionIR.rffSignal;
        % slice selection
        pulseSequence.gzSignal(tgtIdx)      = pulseSequence.gzSignal(tgtIdx) ...
            + repetitionIR.ssSignal;
        % software crushers
        pulseSequence.swcSignal(tgtIdx)     = repetitionIR.swcSignal;
        % time: correct with current time
        pulseSequence.time(tgtIdx)          = repetitionIR.time + currentTime;
        % parts
        pulseSequence.partType(tgtPart)     = repetitionIR.partType;
        pulseSequence.partLimits(tgtPart,:) = repetitionIR.partLimits + idxShift;
        % increase shifts
        idxShift  = idxShift  + repetitionIR.numSteps;
        partShift = partShift + repetitionIR.numParts; 
        % new current time : correct with waiting time for TI
        currentTime = currentTime + irTime;    
        % correct the last time point to be the end time
        pulseSequence.time(idxShift) = currentTime;
    else
        % No IR
        if rr > 1
            targetStartTime = triggerTime + currentCycle*cardiacShotTR;
            % set as time for the current index (last point of prev bssfp)
            % a time just before starting the next bSSFP:
            % avoid potential large RF steps
            pulseSequence.time(idxShift) = targetStartTime - 1e-8;
        end
    end
    % mark the trigger time to sync all bSSFP
    if rr == 1
        triggerTime = currentTime;
    end
    % add all the bSSFP shots of the scheme
    for ss = 1:molliScheme{rr}.shots
        % current shot and verify change of polarity
        currentShot = currentShot + 1;
        actualTI(1,currentShot) = molliScheme{rr}.TI + ...
            (ss-1)*round(cardiacCycleDurationMsec/1e3);
        if rem(currentShot+1,2) && ~evenBssfpReps
            rfScalingFactor = -1;
        else
            rfScalingFactor =  1;
        end
        % check correct timings
        targetStartTime = triggerTime + currentCycle*cardiacShotTR;
        if endShotTime > targetStartTime
            % Throw error: we are computing things badly
            msg = sprintf( ['Time end of prev. bssfp  (%.2fms) larger',...
            ' than target time (%.2fms) for cycle %d'],...
                endShotTime*1e3, targetStartTime*1e3, currentCycle );
            ME = MException('error:edutool', '%s', msg);
            throw(ME);
        end
        % starting time multiple of cardiac cycle
        currentTime = targetStartTime;
         % target indexes
        tgtIdx  = (1:bssfpSequence.numSteps) + idxShift;
        tgtPart = (1:bssfpSequence.numParts) + partShift;
        tgtLimit= (1:bssfpSequence.numEnc)   + limitShift;
        % assign waveforms
        % gradients
        pulseSequence.gxSignal(tgtIdx)      = bssfpSequence.gxSignal;
        pulseSequence.gySignal(tgtIdx)      = bssfpSequence.gySignal;
        pulseSequence.gzSignal(tgtIdx)      = bssfpSequence.gzSignal;
        % rf
        pulseSequence.rfmSignal(tgtIdx)     = bssfpSequence.rfmSignal * rfScalingFactor;
        pulseSequence.rfpSignal(tgtIdx)     = bssfpSequence.rfpSignal;
        pulseSequence.rffSignal(tgtIdx)     = bssfpSequence.rffSignal;
        % readouts and swc
        pulseSequence.swcSignal(tgtIdx)     = bssfpSequence.swcSignal;
        pulseSequence.rxSignal(tgtIdx)      = bssfpSequence.rxSignal;
        % time: correct with current time
        pulseSequence.time(tgtIdx)          = bssfpSequence.time + currentTime;
        % parts
        pulseSequence.partType(tgtPart)     = bssfpSequence.partType;
        pulseSequence.partLimits(tgtPart,:) = bssfpSequence.partLimits + idxShift;
        % rx limits
        pulseSequence.rxLimits(tgtLimit,:)  = bssfpSequence.rxLimits + idxShift;
        % increase shifts
        idxShift  = idxShift  + bssfpSequence.numSteps + 1; % extra step
        partShift = partShift + bssfpSequence.numParts; 
        limitShift= limitShift+ bssfpSequence.numEnc;
        % store last time of the current bssfp to verify correctness
        endShotTime = currentTime + bssfpSequence.time(end);
        % set as time for the current index (last point of this bssfp)
        % a time just before starting the next bSSFP:
        % avoid potential large RF steps
        currentTime = triggerTime + (currentCycle+1)*cardiacShotTR;
        currentTime = currentTime - 1e-8;
        pulseSequence.time(idxShift) = currentTime;
        % correct last part limit to account for extra point
        pulseSequence.partLimits(partShift,2) = pulseSequence.partLimits(partShift,2) + 1;
        % increase cycle
        currentCycle = currentCycle + 1;
    end
    % increas the cycles by the pause
    currentCycle = currentCycle + molliScheme{rr}.pause;
    
end
    
%% add remaining data
pulseSequence.rxSignal(pulseSequence.rxSignal>0) = 1:nnz(pulseSequence.rxSignal);
pulseSequence.timeDiff = reshape(diff([0; pulseSequence.time(:)]),[],1);

%% pulse sequences info
pulseSequence.name          = 'MOLLI';
pulseSequence.type          = 'IR-bSSFP';
pulseSequence.endEvent      = 'none'; % indicates what happens at the end
pulseSequence.totalTime     = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps      = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs        = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps       = bssfpSequence.numReps+1;
pulseSequence.numEnc        = bssfpSequence.numEnc; % number of encoding repetitions
pulseSequence.numShots      = bssfpSequence.numReps; % number of shots
pulseSequence.numTL         = 1; % number of reps in shot
pulseSequence.TE            = bssfpSequence.TE; % repetition TE;
pulseSequence.TR            = bssfpSequence.TR; % repetition TR;
pulseSequence.TI            = molliScheme{1}.TI; % preparation TI;
pulseSequence.ESP           =  0; % echo spacing
pulseSequence.effTE         = bssfpSequence.effTE; % effective echo time
pulseSequence.effTR         = bssfpSequence.TR; % effective TR
pulseSequence.effTI         = molliScheme{1}.TI; % effective TI (to center of K space)
% maximum interleaving possible (slices in a TR)
pulseSequence.maxIL         = 1; % BSSFP sequences do not allow interleaving

pulseSequence.MOLLI.schemeMOLLI = schemeMOLLI;
pulseSequence.MOLLI.pauseMOLLI  = pauseMOLLI;
pulseSequence.MOLLI.TIs         = actualTI;
pulseSequence.MOLLI.textScheme  = acquisition.data.MOLLI.scheme;

%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Number steps        %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Number readouts     %d', pulseSequence.numRxs);
    fprintf(fid, '\n  Number encodings    %d', pulseSequence.numEnc);
    fprintf(fid, '\n  Number repetitions  %d', pulseSequence.numReps);
    fprintf(fid, '\n  TI                  %.3fms', 1e3*pulseSequence.TI);
    fprintf(fid, '\n  TE                  %.3fms', 1e3*pulseSequence.TE);
    fprintf(fid, '\n  TR                  %.3fms', 1e3*pulseSequence.TR);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
end


%% plot sequence
% figure();
% plot(pulseSequence.time,pulseSequence.rfmSignal*1e3);
% xlabel('time (s)');
% hold on
% plot(pulseSequence.time,pulseSequence.gxSignal);
% plot(pulseSequence.time,pulseSequence.gySignal);
% plot(pulseSequence.time,pulseSequence.gzSignal);
% plot(pulseSequence.time(pulseSequence.rxSignal>0), ...
%     pulseSequence.gxSignal(pulseSequence.rxSignal>0), 'o');
% plot(pulseSequence.time(pulseSequence.rxLimits(:,1)), ...
%     pulseSequence.gxSignal(pulseSequence.rxLimits(:,1)), '>');
% plot(pulseSequence.time(pulseSequence.rxLimits(:,2)), ...
%     pulseSequence.gxSignal(pulseSequence.rxLimits(:,2)), '<');
% plot(pulseSequence.time(pulseSequence.partLimits(:,1)),zeros(pulseSequence.numParts,1), '^');
% plot(pulseSequence.time(pulseSequence.partLimits(:,2)),zeros(pulseSequence.numParts,1), 'v');
% if nnz(pulseSequence.swcSignal) > 0
%     plot(pulseSequence.time(pulseSequence.swcSignal>0), ...
%         zeros(nnz(pulseSequence.swcSignal),1), 's', 'LineWidth', 2, 'MarkerSize', 10);
%     legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end', 'SWC');
% else
%     legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end');
% end
