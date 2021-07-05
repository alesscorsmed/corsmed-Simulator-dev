function [pulseSequence] = spinEchoEPI(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYEPI.SPINECHOEPI
%
%	Generates SE EPI sequence.
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
functionName = 'sequence:familyEPI:spinEchoEPI';

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
mainRF      = acquisition.mainRF;
refRF       = acquisition.refRF;
encPlan     = encoding.plan;

%% generate EPI encoding
[timeEPI,signalFE,signalPE,signalRX,~,effTeEPI,effESP,rxLimits,numENC] = ...
    sequence.encodings.cartesianEPI( ...
    acqData.rxBW, acqData.fovFE, acqData.fovPE, ...
    acqData.numFE, acqData.numPE, ...
    acqData.samplingFactorFE, acqData.samplingFactorPE,...
    encPlan.fFactorFE, encPlan.fFactorPE, encPlan.rFactorPE, ...
    mrSystem.maxGStrenght, mrSystem.SlewRate, acqData.gamma, ...
    expControl.sequence.dtGR, expControl);

%% generate the waveforms for the RF 90 (use main RF info)
[rf90Time,rfm90,rfp90,rff90,grSlice90,rfLimits90] = ...
    sequence.excitations.sliceSelectRF( mainRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl );

%% generate the waveforms for the RF 180 (use refocusing RF info)
[rf180Time,rfm180,rfp180,rff180,grSlice180,rfLimits180] = ...
    sequence.excitations.sliceSelectRF( refRF,...
    mrSystem.maxGStrenght, mrSystem.SlewRate, ...
    acqData.gamma, expControl.sequence.dtRF, expControl);

%% verify feasibility of TE and correctness of timing
minTE90  = rf90Time(end)/2 + rf180Time(end)/2;
tWait90  = acqData.TE - minTE90;
minTE180 = rf180Time(end)/2 + effTeEPI;
tWait180 = acqData.TE - minTE180;
minTE = max(minTE180,minTE90);
if tWait90 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        acqData.TE*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
else
    tWait90 = max(tWait90,1e-12);
end
if tWait180 < 0
    msg = sprintf( ['The selected Echo Time (%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is %.2fms'],...
        acqData.TE*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
else
    tWait180 = max(tWait180,1e-12);
end
      
%% combine the into a sequence
numTimesRF90    = length(rf90Time);
numTimesRF180   = length(rf180Time);
numTimesEPI     = length(timeEPI);


% start/end of RF90
iniRF90         = 1;
endRF90         = iniRF90 + numTimesRF90-1;
rfLimits90      = rfLimits90 + iniRF90-1;
% single point after RF90
idx1Wait90      = endRF90 + 1;
% start/end of RF180
iniRF180        = idx1Wait90 + 1;
endRF180        = iniRF180 + numTimesRF180-1;
rfLimits180     = rfLimits180 + iniRF180-1;
% single point after RF180
idxWait180      = endRF180 + 1;
% start/end of EPI encoding
iniEPI          = idxWait180 + 1;
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
pulseSequence.gdwSignal = []; % diffusion gradients

%% assign blocks
pulseSequence.time(iniRF90:endRF90)     = rf90Time;
pulseSequence.time(idx1Wait90)          = tWait90 + pulseSequence.time(idx1Wait90-1);
pulseSequence.time(iniRF180:endRF180)   = rf180Time + pulseSequence.time(iniRF180-1);
pulseSequence.time(idxWait180)          = tWait180 + pulseSequence.time(idxWait180-1);
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

%% EPI
pulseSequence.gxSignal(iniEPI:endEPI) = ...
    pulseSequence.gxSignal(iniEPI:endEPI) + signalFE;
pulseSequence.gySignal(iniEPI:endEPI) = ...
    pulseSequence.gySignal(iniEPI:endEPI) + signalPE;
pulseSequence.rxSignal(iniEPI:endEPI) = signalRX;

%% pulse sequences info
pulseSequence.name    	= 'SE-EPI';
pulseSequence.type     	= 'SE-EPI';
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
pulseSequence.effTE     = acqData.TE; % effective echo time
pulseSequence.effTR     = pulseSequence.totalTime; % effective TR
pulseSequence.effTI     =  0; % effective TI (to center of K space)
pulseSequence.rxLimits  = rxLimits; % starting / end indexes of Readouts
% maximum interleaving possible (slices in a TR)
% TODO: fix this when we have actual TR
pulseSequence.maxIL     = floor(pulseSequence.effTR/pulseSequence.TR);

%% for splitting into parts (RF vs Non-RF)
% multiple parts, depending on its characteristics
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
% rest may have readouts
numParts = numParts + 1;
partLimits = [partLimits; rfLimits180(2)+1, pulseSequence.numSteps];
partType{numParts} = 'RO';
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
    fprintf(fid, '\n  Eff. Time Echo      %.3fms', pulseSequence.effTE*1e3);
    fprintf(fid, '\n  Eff. Echo Spacing   %.3fms', pulseSequence.ESP*1e3 );
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


