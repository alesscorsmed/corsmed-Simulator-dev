function pulse_sequence_IR = generateIRpulse(structIR,structExper)

if strcmp(structIR.type,'sinc')
    if ~isfield(structIR,'IRcycles')
        structIR.IRcycles = 2;
    end
    [B1,~,~] = ...
        pulseSequenceGenerator.addRFsinc(structIR.IRduration,...
        structIR.IRcycles,structIR.angle,...
        structExper.dt,structExper.gamma);
    
    Dph = zeros(size(B1));
    Dfr = zeros(size(B1));    
    
elseif strcmp(structIR.type,'HS1')
    B1max   = 0.0000270606;
    Dwmax   = 1343.39;
    
    [B1,Dph,Dfr] = pulseSequenceGenerator.addRF_AdiabaticHypSec(B1max,...
        Dwmax,structIR.IRduration,structExper.dt);
    
elseif strcmp(structIR.type,'TanTanh')
    
    B1max   = 0.000015;
    fmax    = 9500;
    z       = 10;
    k       = 1.525373047373320;
    [B1,Dph,Dfr] = pulseSequenceGenerator.addRF_AdiabaticTanTanh_FullPassage(B1max,...
        fmax,z,k,structIR.IRduration,structExper.dt);
    
elseif strcmp(structIR.type,'BIR4')
    
    % B1max and Dwmax for a 5ms IR pulse based on the files AM_BIR4.txt and
    % FM_BIR4.txt at your personal documents. It doesn't work for other IR
    % durations.
    B1max   = 0.000015;
    Dwmax   = 8489;
    beta    = 5;
    kappa   = atan(10); %1.4711;
    
    [B1,Dph,Dfr] = pulseSequenceGenerator.addRF_AdiabaticBIR4(B1max,...
        Dwmax,beta,kappa,structIR.angle,structIR.IRduration,structExper.dt);
    
end

if isfield(structIR,'phase')
    Dph = Dph + (structIR.phase)*pi/180; % degrees to rads
end

pulse_sequence_IR.pulse_sequence   = [B1;Dph;Dfr;zeros(3,size(B1,2));repmat(1:size(B1,2),2,1)];
pulse_sequence_IR.isInKspace       = zeros(1,size(B1,2));
pulse_sequence_IR.soft_crushers    = zeros(1,size(B1,2));
pulse_sequence_IR.N_pulse          = size(B1,2);
pulse_sequence_IR.kspace_times     = [];
pulse_sequence_IR.TR_times         = [];
pulse_sequence_IR.receiver_phase   = [];