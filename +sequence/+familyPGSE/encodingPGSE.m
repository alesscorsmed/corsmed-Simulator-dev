function [pulseSequence] = encodingPGSE( ...
    mainRF, refRF, encPG, mrSystem, expControl) 
%
% SEQUENCE.PGSE.ENCODINGPGSE
%
%	Generates PGSE encoding, with 90 and 180 RFs.
%
% INPUT
%
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.pgse.encodingPGSE';

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

gamma = 42.577478518e6;
dt = expControl.sequence.dtRF;

%% generate Pulsed Gradient
[pgTime,pgSignal,pgArea,~,pgPlatLimits] = sequence.waveforms.grTrapPlateau(...
        encPG.TG, encPG.AG, mrSystem.SlewRate, dt);
% compute remaining data    
encPG.Area = pgArea; % gradient area
encPG.Time = pgTime(end); % total time
encPG.TP   = diff(pgTime(pgPlatLimits)); % plateau time
encPG.Beta = (2*pi*gamma)^2 *(encPG.Area^2) *(encPG.TAU-encPG.TG/3);

%% generate the waveforms for the RF 90 (use main RF info)
[rf90Time,rfm90,rfp90,rff90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    gamma, expControl.sequence.dtRF, expControl );

%% generate the waveforms for the RF 180 (use refocusing RF info)
[rf180Time,rfm180,rfp180,rff180,grSlice180,rfLimits180] = ...
    sequence.excitations.sliceSelectRF( refRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    gamma, expControl.sequence.dtRF, expControl );

%% verify feasibility of TE and correctness of timing
minTE90  = rf90Time(end)/2 + pgTime(end) + rf180Time(end)/2;
tWait90  = encPG.TAU - minTE90;
minTE180 = rf180Time(end)/2+ pgTime(end);
tWait180 = encPG.TAU - minTE180;
minTE = max(minTE180,minTE90);
if tWait90 < 0
    msg = sprintf( ['The selected Echo Time (TAU=%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is TAU=%.2fms'],...
        encPG.TAU*1e3, minTE*1e3 );
    eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
end
if tWait180 < 0
    msg = sprintf( ['The selected Echo Time (TAU=%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is TAU=%.2fms'],...
        encPG.TAU*1e3, minTE*1e3 );
    eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
end
      
%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRF180   = length(rf180Time);
numTimesGrad    = length(pgTime);
numTimesWait90  = round(tWait90/dt);
numTimesWait180 = round(tWait180/dt);
numTimesTAU     = round(encPG.TAU/dt);

% start/end of RF90
iniRF90         = 1;
endRF90         = iniRF90 + numTimesRF90-1;
rfLimits90      = rfLimits90 + iniRF90-1;
% start/end of 1st gradient
iniFwGr         = endRF90 + 1;
endFwGr         = iniFwGr+numTimesGrad-1;
% start/end of RF180
iniRF180        = endFwGr + numTimesWait90 + 1;
endRF180        = iniRF180 + numTimesRF180-1;
rfLimits180     = rfLimits180 + iniRF180-1;
% start/end of rewind gradient
iniRwGr         = endRF180 + 1;
endRwGr         = iniRwGr+numTimesGrad-1;
% total time of sequence
numSteps        = max( round(iniRF90 + numTimesRF90/2 + 2*numTimesTAU), ...
    endRF180 + numTimesWait180 );


%% pulse sequence
[pulseSequence] = data.simulation.initializeSequence();

% waveforms
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.timeDiff  = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1); % receiver readout
pulseSequence.swcSignal = zeros(numSteps,1); % software crusher
pulseSequence.gxSignal  = zeros(numSteps,1); % x gradient
pulseSequence.gySignal  = zeros(numSteps,1); % y gradient
pulseSequence.gzSignal  = zeros(numSteps,1); % z gradient
pulseSequence.rfmSignal = zeros(numSteps,1); % RF magnitude
pulseSequence.rfpSignal = zeros(numSteps,1); % RF phase
pulseSequence.rffSignal = zeros(numSteps,1); % RF frequency
pulseSequence.gdwSignal = zeros(numSteps,3); % diffusion gradients

%% assign blocks
pulseSequence.time(:)     = dt*(1:numSteps);
pulseSequence.timeDiff(:) = dt;
%% RF
pulseSequence.rfmSignal(iniRF90:endRF90)    = rfm90;
pulseSequence.rfpSignal(iniRF90:endRF90)    = rfp90;
pulseSequence.gzSignal(iniRF90:endRF90)     = grSlice90;
%% Refocusing
pulseSequence.rfmSignal(iniRF180:endRF180)  = rfm180;
pulseSequence.rfpSignal(iniRF180:endRF180)  = rfp180;
pulseSequence.gzSignal(iniRF180:endRF180)   = grSlice180;
%% PG gradients (Fw and Rw)
pulseSequence.encPG = encPG; % info for diffusion
switch lower(pulseSequence.encPG.Dir)
    case 'x'
        pulseSequence.gxSignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gxSignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gxSignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gxSignal(iniRwGr:endRwGr,1) + pgSignal(:);
    case 'y'
        pulseSequence.gySignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gySignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gySignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gySignal(iniRwGr:endRwGr,1) + pgSignal(:);
    case 'z'
        pulseSequence.gzSignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gzSignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gzSignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gzSignal(iniRwGr:endRwGr,1) + pgSignal(:);
    otherwise
        ME = MException('sequence:wrongGradDir',...
            '%s : gradient direction %s not supported',functionName,gradDir);
        throw(ME);
end
pulseSequence.rxSignal(end)    = 1; % single RO at end
rxLimits = [numSteps, numSteps];

%% pulse sequences info
pulseSequence.name         = 'Pulsed Gradient Spin Echo';
pulseSequence.type         = 'PG-SE';
pulseSequence.endEvent     = 'none'; % indicates what happens at the end
pulseSequence.totalTime    = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps     = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs       = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.effESP       = 0;
pulseSequence.effTE        = encPG.TAU;
pulseSequence.effTR        = pulseSequence.totalTime;
pulseSequence.rxLimits     = rxLimits; % starting / end indexes of Readouts

%% for splitting into parts (RF vs Non-RF)
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
% second gradient
numParts = numParts + 1;
partLimits = [partLimits; rfLimits180(2)+1, endRwGr];
partType{numParts} = 'GR';
% rest may have readouts
if pulseSequence.numSteps > endRwGr
    numParts = numParts + 1;
    partLimits = [partLimits; endRwGr+1, pulseSequence.numSteps];
    partType{numParts} = 'RO';
end
% assign
pulseSequence.numParts   = numParts; % number of parts
pulseSequence.partType   = partType; % type of part
pulseSequence.partLimits = partLimits; % index start/end of part

%% update the acquisition data

%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Number steps        %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Number readouts     %d', pulseSequence.numRxs);
    fprintf(fid, '\n  DW Grad. direction  %s', pulseSequence.encPG.Dir );
    fprintf(fid, '\n  DW Grad. Time       %.3fms (Plateau %.3fms)',...
        pulseSequence.encPG.Time*1e3, pulseSequence.encPG.TP*1e3 );
    fprintf(fid, '\n  DW Grad. Amplitude  %.3fmT/m', pulseSequence.encPG.AG*1e3);
    fprintf(fid, '\n  DW Grad. Area       %.3e', pulseSequence.encPG.Area);
    fprintf(fid, '\n  Effective TAU       %.3fms', pulseSequence.encPG.TAU*1e3);
    fprintf(fid, '\n  Effective Beta      %.3f', pulseSequence.encPG.Beta );
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
%         % plot
%         plot(pulseSequence.time,pulseSequence.rfmSignal*1e3)
%         hold on
%         plot(pulseSequence.time,pulseSequence.gxSignal +...
%             pulseSequence.gdwSignal(:,1));
%         plot(pulseSequence.time,pulseSequence.gySignal +...
%             pulseSequence.gdwSignal(:,2));
%         plot(pulseSequence.time,pulseSequence.gzSignal +...
%             pulseSequence.gdwSignal(:,3));

end


