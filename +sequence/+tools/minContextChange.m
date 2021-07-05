function [pulseSequence] = minContextChange(...
    pulseSequence, expControl)
%
% SEQUENCE.TOOLS.MINCONTEXTCHANGE
%
%	Combines sequence parts to 
%   minimize the number of kernels called during simulation.
%
% INPUT
%   pulseSequence    original sequence    
%   expControl      
%
% OUTPUT
%   pulseSequence    new sequence with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.tools.minContextChange';

% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

% allocate for new partition
newNumParts     = 0;
newPartType     = [];
newPartLimits   = [];

% initialize first partition
currentStart = pulseSequence.partLimits(1,1);
currentType  = pulseSequence.partType{1};

% initial number of sequence parts
numParts = pulseSequence.numParts;

% loop on parts and combine if possible
for ss = 2:numParts
    
    %% define the cases when we need to change context
    switch pulseSequence.partType{ss}
        
        case 'RF'
            %% RF as different type
            if strcmpi(currentType, 'rf')
                currentType = 'RF';
                changeContext = 0;
            else
                % For the analytical case:
                %  more aggressive: can combine RFs and GR
                if strcmpi(expControl.simulation.simulationEngine, 'analytical') ...
                        && strcmpi(currentType, 'gr')
                    % combine as RF
                    currentType = 'RF';
                    changeContext = 0;
                else
                    % For Numerial / Phasor / Diffusion, keep separated
                    changeContext = 1;
                end
            end
            
        case 'RO'
            %% can absorb GR
            if strcmpi(currentType, 'ro') ...
                    || strcmpi(currentType, 'gr')
                currentType = 'RO';
                changeContext = 0;
            else
                changeContext = 1;
            end
            
        case 'GR'
            %% pure gradients
            if ~strcmpi(expControl.simulation.simulationEngine, 'analytical') ...
                    && strcmpi(currentType, 'rf')
                    % change context to separate GR from RF 
                    % unless is analytical
                    currentType = 'RF';
                    changeContext = 1;
            else
                % For the analytical case or current type not RF:
                %  GR can be combined with all remaining cases
                changeContext = 0;
                % do not assign currentType, since GR part will be absorbed
            end
            
        otherwise                        
            %% DF case: can absorb GR
            if strcmpi(currentType, 'df') ...
                    || strcmpi(currentType, 'gr')
                currentType = 'DF';
                changeContext = 0;
            else
                changeContext = 1;
            end
    end
    
    %% if we need to change context
    if changeContext
        % create a new part
        newNumParts              = newNumParts + 1;
        newPartType{newNumParts} = currentType;
        newPartLimits            = [newPartLimits; ...
            currentStart, pulseSequence.partLimits(ss-1,2)];
        % update current type
        currentStart = pulseSequence.partLimits(ss,1);
        currentType  = pulseSequence.partType{ss};
    end
    
end

%% add last part
newNumParts              = newNumParts + 1;
newPartType{newNumParts} = currentType;
newPartLimits            = [newPartLimits; ...
    currentStart, pulseSequence.partLimits(ss,2)];

%% assign to sequence
pulseSequence.numParts   = newNumParts;
pulseSequence.partType   = newPartType;   
pulseSequence.partLimits = newPartLimits;


%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  Original   # parts  %d', numParts);
    fprintf(fid, '\n  Compressed # parts  %d', newNumParts);
    fprintf(fid, '\n  Reduction           %.1f%%',...
        100*(1- newNumParts/numParts));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

