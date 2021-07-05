function pulseSequence = convertOldToNewPulseSequence(pulseSequence,...
    pulseSeqFamilyName,pulse_sequence,isInKspace,soft_crushers,dt)

% basic info
pulseSequence.name          = pulseSeqFamilyName;
pulseSequence.type          = pulseSeqFamilyName;
pulseSequence.endEvent      = 'none'; % indicates what happens at the end
pulseSequence.numSteps      = size(pulse_sequence,2);

pulseSequence.time          = zeros(pulseSequence.numSteps,1); % times
pulseSequence.timeDiff      = zeros(pulseSequence.numSteps,1); % time deltas
pulseSequence.rxSignal      = zeros(pulseSequence.numSteps,1); % receiver readout
pulseSequence.swcSignal     = zeros(pulseSequence.numSteps,1); % software crusher
pulseSequence.gxSignal      = zeros(pulseSequence.numSteps,1); % x gradient
pulseSequence.gySignal      = zeros(pulseSequence.numSteps,1); % y gradient
pulseSequence.gzSignal      = zeros(pulseSequence.numSteps,1); % z gradient
pulseSequence.rfmSignal     = zeros(pulseSequence.numSteps,1); % RF magnitude
pulseSequence.rfpSignal     = zeros(pulseSequence.numSteps,1); % RF phase
pulseSequence.rffSignal     = zeros(pulseSequence.numSteps,1); % RF frequency

pulseSequence.timeDiff(:)   = (pulse_sequence(8,:)-...
    pulse_sequence(7,:)+1)*dt;
pulseSequence.time(:)       = cumsum(pulseSequence.timeDiff);

pulseSequence.rxSignal(:)   = isInKspace(:);
pulseSequence.swcSignal(:)  = soft_crushers(:);

pulseSequence.rfmSignal(:)  = pulse_sequence(1,:);
pulseSequence.rfpSignal(:)	= pulse_sequence(2,:);
pulseSequence.rffSignal(:)  = pulse_sequence(3,:);
pulseSequence.gxSignal(:)   = pulse_sequence(4,:);
pulseSequence.gySignal(:)	= pulse_sequence(5,:);
pulseSequence.gzSignal(:)   = pulse_sequence(6,:);

pulseSequence.totalTime     = pulseSequence.time(end);
pulseSequence.numRxs        = nnz(pulseSequence.rxSignal(:));

pulseSequence.numParts      = 1; % number of parts
pulseSequence.partType{1}   = 'DF'; % type of part: RF / RO / GR / DF
pulseSequence.partLimits    = [1, pulseSequence.numSteps]; % index start/end of part
pulseSequence.rxLimits      = [1, pulseSequence.numSteps];
