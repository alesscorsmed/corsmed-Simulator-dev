function [pulseSequence] = pulsedGradientSpinEchoEPI(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYPGSE.PULSEDGRADSPINECHOEPI
%
%	Generates PGSE EPI sequence.
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
functionName = 'sequence:familyPGSE:pulsedGradientSpinEchoEPI';

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
encPG       = acquisition.encPG;
mainRF      = acquisition.mainRF;
refRF       = acquisition.refRF;
encPlan     = encoding.plan;

%% generate EPI encoding
[timeEPI,signalFE,signalPE,signalRX,dtEPI,effTeEPI,effESP,rxLimits,numENC] = ...
    sequence.encodings.cartesianEPI( ...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, ...
    acqData.numFE, acqData.numPE, ...
    acqData.samplingFactorFE, acqData.samplingFactorPE,...
    encPlan.fFactorFE, encPlan.fFactorPE, encPlan.rFactorPE, ...
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma, ...
    expControl.sequence.dtGR, expControl );

%% generate Pulsed Gradient
[pgTime,pgSignal,pgArea,~,pgPlatLimits] = sequence.waveforms.grTrapPlateau(...
        encPG.TG, encPG.AG, mrSystem.SlewRate, dtEPI);
% compute remaining data    
encPG.Area = pgArea; % gradient area
encPG.Time = pgTime(end); % total time
encPG.TP   = diff(pgTime(pgPlatLimits)); % plateau time
encPG.Beta = (2*pi*acqData.gamma)^2 *(encPG.Area^2) *(encPG.TAU-encPG.TG/3);

%% generate the waveforms for the RF 90 (use main RF info)
[rf90Time,rfm90,rfp90,rff90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% generate the waveforms for the RF 180 (use refocusing RF info)
[rf180Time,rfm180,rfp180,rff180,grSlice180,rfLimits180] = ...
    sequence.excitations.sliceSelectRF( refRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% verify feasibility of TE and correctness of timing
minTE90  = rf90Time(end)/2 + pgTime(end) + rf180Time(end)/2;
tWait90  = encPG.TAU - minTE90;
minTE180 = rf180Time(end)/2+ pgTime(end) + effTeEPI;
tWait180 = encPG.TAU - minTE180;
minTE = max(minTE180,minTE90);
if tWait90 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        encPG.TAU*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end
if tWait180 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        encPG.TAU*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end
      
%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRF180   = length(rf180Time);
numTimesGrad    = length(pgTime);
numTimesEPI     = length(timeEPI);
numTimesWait90  = round(tWait90/dtEPI);
numTimesWait180 = round(tWait180/dtEPI);

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
% start/end of EPI encoding
iniEPI          = endRwGr + numTimesWait180 + 1;
endEPI          = iniEPI+numTimesEPI-1;
% shift rx Limits of the encoding
rxLimits        = rxLimits + iniEPI-1;
% total time of sequence
numSteps        = endEPI;

% to define the slice selection
if expControl.sequence.deactivateSS
    ssScale = 0; % zero slice selection grad
else
    ssScale = 1;
end

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
pulseSequence.time(iniRF90:endRF90)     = rf90Time;
pulseSequence.time(iniFwGr:endFwGr)     = pgTime  + pulseSequence.time(iniFwGr-1);
pulseSequence.time(endFwGr+1:iniRF180-1)= tWait90 + pulseSequence.time(endFwGr);
pulseSequence.time(iniRF180:endRF180)   = rf180Time + pulseSequence.time(iniRF180-1);
pulseSequence.time(iniRwGr:endRwGr)     = pgTime + pulseSequence.time(iniRwGr-1);
pulseSequence.time(endRwGr+1:iniEPI-1)  = tWait180 + pulseSequence.time(endRwGr);
pulseSequence.time(iniEPI:endEPI)       = timeEPI + pulseSequence.time(iniEPI-1);
pulseSequence.timeDiff(:)               = diff([0; pulseSequence.time]);
%% RF
pulseSequence.rfmSignal(iniRF90:endRF90)    = rfm90;
pulseSequence.rfpSignal(iniRF90:endRF90)    = rfp90;
pulseSequence.rffSignal(iniRF90:endRF90)    = rff90;
pulseSequence.gzSignal(iniRF90:endRF90)     = ssScale*grSlice90;
%% Refocusing
pulseSequence.rfmSignal(iniRF180:endRF180)  = rfm180;
pulseSequence.rfpSignal(iniRF180:endRF180)  = rfp180;
pulseSequence.rffSignal(iniRF180:endRF180)  = rff180;
pulseSequence.gzSignal(iniRF180:endRF180)   = ssScale*grSlice180;
%% PG gradients (Fw and Rw)
pulseSequence.encPG = encPG; % info for diffusion
switch lower(pulseSequence.encPG.Dir)
    case 'x'
        pulseSequence.gdwSignal(iniFwGr:endFwGr,1) =...
            pulseSequence.gdwSignal(iniFwGr:endFwGr,1) + pgSignal(:);
        pulseSequence.gdwSignal(iniRwGr:endRwGr,1) =...
            pulseSequence.gdwSignal(iniRwGr:endRwGr,1) + pgSignal(:);
    case 'y'
        pulseSequence.gdwSignal(iniFwGr:endFwGr,2) =...
            pulseSequence.gdwSignal(iniFwGr:endFwGr,2) + pgSignal(:);
        pulseSequence.gdwSignal(iniRwGr:endRwGr,2) =...
            pulseSequence.gdwSignal(iniRwGr:endRwGr,2) + pgSignal(:);
    case 'z'
        pulseSequence.gdwSignal(iniFwGr:endFwGr,3) =...
            pulseSequence.gdwSignal(iniFwGr:endFwGr,3) + pgSignal(:);
        pulseSequence.gdwSignal(iniRwGr:endRwGr,3) =...
            pulseSequence.gdwSignal(iniRwGr:endRwGr,3) + pgSignal(:);
    otherwise
        ME = MException('sequence:wrongGradDir',...
            '%s : gradient direction %s not supported',functionName,gradDir);
        throw(ME);
end
%% EPI
pulseSequence.gxSignal(iniEPI:endEPI) = ...
    pulseSequence.gxSignal(iniEPI:endEPI) + signalFE;
pulseSequence.gySignal(iniEPI:endEPI) = ...
    pulseSequence.gySignal(iniEPI:endEPI) + signalPE;
pulseSequence.rxSignal(iniEPI:endEPI) = signalRX;

%% pulse sequences info
pulseSequence.name    	= 'Pulsed Gradient SE-EPI';
pulseSequence.type    	= 'PG-SE-EPI';
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.totalTime = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps  = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs    = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps   = 1;
pulseSequence.numEnc    = numENC; % number of encoding repetitions
pulseSequence.numShots  = 1; % number of shots
pulseSequence.numTL     = numENC; % number of reps in shot
pulseSequence.TE        = effTeEPI; % repetition TE;
pulseSequence.TR        = pulseSequence.totalTime; % repetition TR;
pulseSequence.TI        =  0; % preparation TI;
pulseSequence.ESP       = effESP; % echo spacing
pulseSequence.effTE     = encPG.TAU; % effective echo time
pulseSequence.effTR     = pulseSequence.totalTime; % effective TR
pulseSequence.effTI     =  0; % effective TI (to center of K space)
pulseSequence.rxLimits 	= rxLimits; % starting / end indexes of Readouts
% maximum interleaving possible (slices in a TR)
% TODO: fix this when we have actual TR
pulseSequence.maxIL     = floor(pulseSequence.effTR/pulseSequence.TR);

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
    
    %     % plot
    %     plot(pulseSequence.time,pulseSequence.rfmSignal*1e3)
    %     hold on
    %     plot(pulseSequence.time(rfLimits90),[0,0], '^');
    %     plot(pulseSequence.time(rfLimits180),[0,0], 'v');
    %     plot(pulseSequence.time,pulseSequence.gxSignal);
    %     plot(pulseSequence.time,pulseSequence.gySignal);

end


