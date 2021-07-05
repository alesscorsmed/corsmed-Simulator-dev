function [pulse_sequence,N_pulseIR,N_IRmodule] = IRgenerator(dt,gamma,...
    TI,IRduration,IRtype,TIendPoint,conn_localdb,...
    experiment_id,pulseq_id,pulseSeqFamilyName)
% N_pulseIR holds the length of the pulse_sequence, after the application
% of the fast algorithm. The IRgenerator generates the IR module with
% fast==1 as default option.
% N_TI holds the length of the pulse sequence in real life

IRcycles    = 2;
IRangle     = 180;

% STEP 1 - Design RF
if strcmp(IRtype,'sinc')    
    [RF,BW,rf_timesteps] = pulseSequenceGenerator.addRFsinc(IRduration,IRcycles,...
        IRangle,dt,gamma); %#ok
end

if strcmp(pulseSeqFamilyName,'IR-bSSFP') || strcmp(pulseSeqFamilyName,'bSSFP')
    % In bSSFP the inversion time counts from the end of the IR pulse till 
    % the middle of the central kspace line of the bSSFP acquisition. N_TI 
    % holds the total timesteps from the end of the IR pulse till the start
    % of the bSSFP readout.
    N_TI = floor(TI/dt) - TIendPoint;
elseif strcmp(pulseSeqFamilyName,'TSE') || ...
        strcmp(pulseSeqFamilyName,'IR-TSE') || ...
        strcmp(pulseSeqFamilyName,'SS-FSE')
    % In TSE the inversion time counts from the mid of the IR pulse till 
    % the mid of the 90 degree RF pulse.
    N_TI = floor(TI/dt) - TIendPoint - round(size(RF,2)/2);
end

if N_TI<=0
    if strcmp(pulseSeqFamilyName,'IR-bSSFP') || strcmp(pulseSeqFamilyName,'bSSFP')
        msg = ['The selected Inversion Time (TI) is smaller than the minimum available TI (',...
        num2str(TIendPoint*dt),'sec) for the current configuration. Consider',...
        ' increasing the TI or decreasing the duration of the bSSFP readout.'];
    elseif strcmp(pulseSeqFamilyName,'TSE') || ...
            strcmp(pulseSeqFamilyName,'IR-TSE') || ...
            strcmp(pulseSeqFamilyName,'SS-FSE')
        msg = ['The selected Inversion Time (TI) is smaller than the minimum available TI (',...
        num2str((TIendPoint+round(size(RF,2)/2))*dt),'sec) for the current configuration. Consider',...
        ' increasing the TI, decreasing the duration of the IR pulse or the duration',...
        ' of the 90 degree RF pulse.'];
    end
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end

% For bSSFP, N_IRmodule holds the total timesteps of the IR module (IR pulse +
% remaining till the start of the bSSFP readout). For TSE, N_IRmodule holds
% the total timesteps of the IR module + remaining timesteps till the start
% of the 90 degree RF pulse.
N_IRmodule = rf_timesteps + N_TI;

% N_pulseIR holds the length of the compressed version (fast algorithm)...
% of the IR module
N_pulseIR  = rf_timesteps + 1;

pulse_sequence                      = zeros(8,N_pulseIR);
pulse_sequence(1,1:rf_timesteps)    = RF;
pulse_sequence(7,:)                 = 1:size(pulse_sequence,2);
pulse_sequence(8,:)                 = 1:size(pulse_sequence,2);

% Add one extra point that holds the accumulated effect of the pulse 
% sequence till the end of the IR module
pulse_sequence(8,end) = N_IRmodule;