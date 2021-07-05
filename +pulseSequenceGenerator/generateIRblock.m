function [pulse_sequence_IRmodule,blockActualTimesteps,blockCompressedTimesteps] = ...
    generateIRblock(structIR,structExper,TIendPoint,pulseSeqFamilyName,...
    compressIRmodule,conn_localdb,experiment_id,pulseq_id)
% INPUTS:
% TIendPoint counts the number of points of the host pulse sequence that
% are considered part of TI.
% OUTPUTS:
% blockActualTimesteps holds the "true" number of timesteps of the
% blockIR. blockCompressedTimesteps holds the "effective" number of 
% timesteps due to the application of the fast algorithm.

pulse_sequence_IR = pulseSequenceGenerator.generateIRpulse(structIR,structExper);

if strcmp(pulseSeqFamilyName,'TSE') || ...
        strcmp(pulseSeqFamilyName,'IR-TSE') || ...
        strcmp(pulseSeqFamilyName,'SS-FSE') || ...
        strcmp(pulseSeqFamilyName,'MP-RAGE')
    % In these pulse sequence types the inversion time counts from the mid 
    % of the IR pulse till the mid of the 1st excitation RF pulse.
    N_TI = floor(structIR.TI/structExper.dt) - TIendPoint - ...
        round(size(pulse_sequence_IR.pulse_sequence,2)/2);
else
    % In bSSFP, the inversion time counts from the end of the IR pulse till 
    % the middle of the central kspace line of the bSSFP acquisition. N_TI 
    % holds the total timesteps from the end of the IR pulse till the start
    % of the bSSFP readout.
    % In SE, the inversion time counts from the end of the IR pulse till
%     the start of the host pulse sequence. In this case, TIendPoint is 0
    N_TI = floor(structIR.TI/structExper.dt) - TIendPoint;
end

if N_TI<=0
    if strcmp(pulseSeqFamilyName,'TSE') || ...
            strcmp(pulseSeqFamilyName,'IR-TSE') || ...
            strcmp(pulseSeqFamilyName,'SS-FSE')
        msg = ['The selected Inversion Time (TI) is smaller than the minimum available TI (',...
        num2str((TIendPoint+round(size(pulse_sequence_IR.pulse_sequence,2)/2))*structExper.dt),'sec) for the current configuration. Consider',...
        ' increasing the TI, decreasing the duration of the IR pulse or the duration',...
        ' of the 90 degree RF pulse.'];
    else
        msg = ['The selected Inversion Time (TI) is smaller than the minimum available TI (',...
        num2str(TIendPoint*structExper.dt),'sec) for the current configuration. Consider',...
        ' increasing the TI or decreasing the duration of the ',...
        pulseSeqFamilyName,' readout.'];
    end
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end

% blockActualTimesteps holds the total timesteps of the IR module. For bSSFP, 
% blockActualTimesteps holds the IR pulse + remaining till the start of the 
% bSSFP readout. For TSE, blockActualTimesteps holds the total timesteps of 
% the IR module + remaining timesteps till the start of the 90 degree RF pulse.
blockActualTimesteps = size(pulse_sequence_IR.pulse_sequence,2) + N_TI;

% blockCompressedTimesteps holds the length of the compressed version (fast 
% algorithm) of the IR module
if TIendPoint == 0
    blockCompressedTimesteps  = size(pulse_sequence_IR.pulse_sequence,2);
else
    blockCompressedTimesteps  = size(pulse_sequence_IR.pulse_sequence,2) + 1;
end

if compressIRmodule
    pulse_sequence_IRmodule         = zeros(8,blockCompressedTimesteps);
    
    pulse_sequence_IRmodule(:,1:size(pulse_sequence_IR.pulse_sequence,2)) = ...
        pulse_sequence_IR.pulse_sequence;
    
    pulse_sequence_IRmodule(7,:)    = 1:size(pulse_sequence_IRmodule,2);
    pulse_sequence_IRmodule(8,:)    = 1:size(pulse_sequence_IRmodule,2);
    
    % Add one extra point that holds the accumulated effect of the pulse 
    % sequence till the end of the IR module
    pulse_sequence_IRmodule(8,end)  = blockActualTimesteps;
else
    pulse_sequence_IRmodule         = zeros(8,blockActualTimesteps);
    
    pulse_sequence_IRmodule(:,1:size(pulse_sequence_IR.pulse_sequence,2)) = ...
        pulse_sequence_IR.pulse_sequence;
    
    pulse_sequence_IRmodule(7,:)    = 1:size(pulse_sequence_IRmodule,2);
    pulse_sequence_IRmodule(8,:)    = 1:size(pulse_sequence_IRmodule,2);
end