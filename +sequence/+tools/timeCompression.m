function [compressed] = timeCompression(...
    repetition, expControl)
%
% SEQUENCE.TOOLS.TIMECOMPRESSION
%
%	Applies time compression of a repetion based on accrued area.
%
% INPUT
%   repetition(iRep)  original repetition(iRep)    
%   expControl      
%
% OUTPUT
%   compressed(iRep)  compressed(iRep) data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.tools.timeCompression';

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

for iRep = size(repetition,2):-1:1
    %% keep entries that are critical: RF, RX, SWC
    idxAll= 1:repetition(iRep).numSteps;
    idxKeep = idxAll( (repetition(iRep).rfEntries > 0) ...
        | (repetition(iRep).rxSignal > 0) | (repetition(iRep).swcSignal > 0) );
    % add limits of the different parts
    idxKeep = union(idxKeep, repetition(iRep).partLimits(:));
    numSteps = length(idxKeep);

    %% new data
    compressed(iRep) = repetition(iRep);
    compressed(iRep).numSteps = numSteps;
    % allocation of waveforms
    compressed(iRep).time      = zeros(numSteps,1);
    compressed(iRep).rxSignal  = zeros(numSteps,1); % receiver readout
    compressed(iRep).swcSignal = zeros(numSteps,1); % software crusher
    compressed(iRep).feSignal  = zeros(numSteps,1); % freq encoding (x) gradient
    compressed(iRep).peSignal  = zeros(numSteps,1); % phase encoding (y) gradient
    compressed(iRep).seSignal  = zeros(numSteps,1); % 3D phase (slice) encoding (z) gradient
    compressed(iRep).rfmSignal = zeros(numSteps,1); % RF magnitude
    compressed(iRep).rfpSignal = zeros(numSteps,1); % RF phase
    compressed(iRep).rffSignal = zeros(numSteps,1); % RF frequency
    compressed(iRep).ssSignal  = zeros(numSteps,1); % slice selection gradient
    compressed(iRep).rfEntries = zeros(numSteps,1); % ones or zeros depending on RF

    %% unchanged waveforms
    compressed(iRep).time(:)      = repetition(iRep).time(idxKeep);
    compressed(iRep).rxSignal(:)  = repetition(iRep).rxSignal(idxKeep);
    compressed(iRep).swcSignal(:) = repetition(iRep).swcSignal(idxKeep);
    compressed(iRep).rfmSignal(:) = repetition(iRep).rfmSignal(idxKeep);
    compressed(iRep).rfpSignal(:) = repetition(iRep).rfpSignal(idxKeep);
    compressed(iRep).rffSignal(:) = repetition(iRep).rffSignal(idxKeep);
    compressed(iRep).rfEntries(:) = repetition(iRep).rfEntries(idxKeep);

    %% compress the gradient waveforms by area
    % cumulative area, first entry unchanged
    tDiff   = [repetition(iRep).time(1); diff(repetition(iRep).time)];
    feArea  = cumsum(repetition(iRep).feSignal.*tDiff);
    peArea  = cumsum(repetition(iRep).peSignal.*tDiff);
    seArea  = cumsum(repetition(iRep).seSignal.*tDiff);
    ssArea  = cumsum(repetition(iRep).ssSignal.*tDiff);
    % new signals are:
    %   area increase for each entry in the compressed(iRep) format
    %   divided by the time increment in the compressed(iRep) format
    %   note first entry unchanged
    % we are keeping larger time steps, with same accrued phase (area)
    tDiff                   = diff(compressed(iRep).time);
    compressed(iRep).feSignal(:)  = [repetition(iRep).feSignal(1); diff(feArea(idxKeep))./tDiff];
    compressed(iRep).peSignal(:)  = [repetition(iRep).peSignal(1); diff(peArea(idxKeep))./tDiff];
    compressed(iRep).seSignal(:)  = [repetition(iRep).seSignal(1); diff(seArea(idxKeep))./tDiff];
    compressed(iRep).ssSignal(:)  = [repetition(iRep).ssSignal(1); diff(ssArea(idxKeep))./tDiff];

    %% find new part limits
    compressed(iRep).numParts     = repetition(iRep).numParts;
    compressed(iRep).partType     = repetition(iRep).partType;
    compressed(iRep).partLimits   = repetition(iRep).partLimits;
    for rr = 1:repetition(iRep).numParts
        compressed(iRep).partLimits(rr,1) = find(idxKeep == repetition(iRep).partLimits(rr,1));
        compressed(iRep).partLimits(rr,2) = find(idxKeep == repetition(iRep).partLimits(rr,2));
    end

    %% find new readout limits
    for rr = 1:size(repetition(iRep).rxLimits,1)
        compressed(iRep).rxLimits(rr,1) = find(idxKeep == repetition(iRep).rxLimits(rr,1));
        compressed(iRep).rxLimits(rr,2) = find(idxKeep == repetition(iRep).rxLimits(rr,2));
    end

    % verify areas keep consistent
    tDiffCmp  = [compressed(iRep).time(1); diff(compressed(iRep).time)];
    tDiffDef  = [repetition(iRep).time(1); diff(repetition(iRep).time)];

    feError   = sum(compressed(iRep).feSignal.*tDiffCmp) ...
        - sum(repetition(iRep).feSignal.*tDiffDef);
    peError   = sum(compressed(iRep).peSignal.*tDiffCmp) ...
        - sum(repetition(iRep).peSignal.*tDiffDef);
    seError   = sum(compressed(iRep).seSignal.*tDiffCmp) ...
        - sum(repetition(iRep).seSignal.*tDiffDef);
    ssError   = sum(compressed(iRep).ssSignal.*tDiffCmp) ...
        - sum(repetition(iRep).ssSignal.*tDiffDef);

    if abs(feError) > 1e-12
        ME = MException('compression:wrongGradArea',...
            '%s : Freq. Encoding gradient area do not match after compression',functionName);
        throw(ME);
    end
    if abs(peError) > 1e-12
        ME = MException('compression:wrongGradArea',...
            '%s : Phase Encoding gradient area do not match after compression',functionName);
        throw(ME);
    end
    if abs(seError) > 1e-12
        ME = MException('compression:wrongGradArea',...
            '%s : Slice Encoding gradient area do not match after compression',functionName);
        throw(ME);
    end
    if abs(ssError) > 1e-12
        ME = MException('compression:wrongGradArea',...
            '%s : Slice Selection gradient area do not match after compression',functionName);
        throw(ME);
    end

    %% final message
    if  expControl.debug.debugMode
        fprintf(fid, ...
            '\n%s : done for %s repetition (#%d), elapsed time %.3fs',...
            functionName, repetition(iRep).type, iRep, toc(tTotal));
        fprintf(fid, '\n  IRL Time            %.3fs', compressed(iRep).totalTime);
        fprintf(fid, '\n  Original   # steps  %d', repetition(iRep).numSteps);
        fprintf(fid, '\n  compressed # steps  %d', compressed(iRep).numSteps);
        fprintf(fid, '\n  Reduction           %.1f%%',...
            100*(1- compressed(iRep).numSteps/repetition(iRep).numSteps));
        fprintf(fid, '\n');
        if fid ~=1
            fclose(fid);
        end
    end

end