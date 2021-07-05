function [pulseSequence] = reversedPG(...
     encPG, mrSystem, expControl)
%
% SEQUENCE.PGSE.REVERSEDPG
%
%	Generates PGSE encoding, with only gradients, one FW, one RW.
%   Same effect as PGSE but without RFs.
%
% INPUT
%
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.pgse.reversedPG';

if (nargin<3)
    ME = MException('sequence:wrongArgument',...
        '%s : invalid arguments',functionName);
    throw(ME);
end

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
    
%% combine the gradients into a sequence
tWait        = 10e-6;
numTimesGrad = length(pgTime);
numTimesWait = round(tWait/dt);
numTimesTAU  = max(round(encPG.TAU/dt), numTimesGrad+numTimesWait);

% start/end of forward gradient
iniFwGr         = numTimesWait;
endFwGr         = iniFwGr+numTimesGrad-1;
% start/end of rewind gradient
iniRwGr         = iniFwGr+numTimesTAU;
endRwGr         = iniRwGr+numTimesGrad-1;
% total time of sequence
numSteps        = ceil(iniRwGr + numTimesGrad/2 + numTimesTAU);

% update and compute remaining PG data
encPG.TAU  = dt*numTimesTAU;
encPG.Area = pgArea; % gradient area
encPG.Time = pgTime(end); % total time
encPG.TP   = diff(pgTime(pgPlatLimits)); % plateau time
encPG.Beta = (2*pi*gamma)^2 *(encPG.Area^2) *(encPG.TAU-encPG.TG/3);

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

% for area compressed sequences
pulseSequence.hArea     = zeros(numSteps,1); % time deltas for areas
pulseSequence.gxArea    = zeros(numSteps,1); % x gradient area delta
pulseSequence.gyArea    = zeros(numSteps,1); % y gradient area delta
pulseSequence.gzArea    = zeros(numSteps,1); % z gradient area delta

%% assign blocks
pulseSequence.time(:)     = dt*(1:numSteps);
pulseSequence.timeDiff(:) = dt;
%% PG gradients (Fw and Rw)
pulseSequence.encPG = encPG; % info for diffusion
switch lower(pulseSequence.encPG.Dir)
    case 'x'
        pulseSequence.gxSignal(iniFwGr:endFwGr,1) =  pgSignal(:);
        pulseSequence.gxSignal(iniRwGr:endRwGr,1) = -pgSignal(:);
    case 'y'
        pulseSequence.gySignal(iniFwGr:endFwGr,1) =  pgSignal(:);
        pulseSequence.gySignal(iniRwGr:endRwGr,1) = -pgSignal(:);
    case 'z'
        pulseSequence.gzSignal(iniFwGr:endFwGr,1) =  pgSignal(:);
        pulseSequence.gzSignal(iniRwGr:endRwGr,1) = -pgSignal(:);
    otherwise
        ME = MException('sequence:wrongGradDir',...
            '%s : gradient direction %s not supported',functionName,gradDir);
        throw(ME);
end
pulseSequence.rxSignal(:) = 1:numSteps;
rxLimits = [1, numSteps];

%% pulse sequences info
pulseSequence.name         = 'Reversed Pulsed Gradient';
pulseSequence.type         = 'reversedPG';
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
partLimits = [1, numSteps];
partType{1} = 'RO';
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


