function [time,rfm,rfp,rff,gr,rfLimits] = sliceSelectRF(rfPulse,...
    maxGStrenght, gradSlewrate, gamma, dt, expControl)
%
% SEQUENCE.EXCITATIONS.SLICESELECTRF
%
%	Generates a RFBlock with optional slice selection,
%   and pre/post gradient rewinds.
%
% INPUT
%   duration        total RF duration, in s
%   cycles          number of cycles of the sinc
%   angle           Flip angle, in degrees
%   tstep           time discretization
%   gamma           gyromagnetic ratio
%
% OUTPUT
%   time            discretized time vector, starts in tstep
%   signal          signal vector
%   BW              RF bandwidth, in Hz
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.excitations.sliceSelectRF';

if (nargin < 1 || isempty(rfPulse))
    ME = MException('sequence:wrongArgument',...
        '%s : invalid rfPulse structure',functionName);
    throw(ME);
end
if (nargin < 2 || isempty(maxGStrenght))
    maxGStrenght=0.030;
end
if (nargin < 3 || isempty(gradSlewrate))
    gradSlewrate=150.0;
end
if (nargin < 4 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end
if (nargin < 5 || isempty(dt))
    dt = 1e-6;
end
if (nargin < 6 || isempty(expControl))
    expControl.connLocalDB = [];
    expControl.application = 'unknown';
    expControl.debug.debugMode = 0;
    expControl.debug.debugFile = '';
end

% info for debugging
if expControl.debug.debugMode
    % open file if possible, otherwise dump to stdout
    try
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

% STEP 1 - Design RF
switch lower(rfPulse.type)
    case 'hs1'
        rfPulse.doSliceSelect = 0;
        % predefined values
        b1Max   = 0.0000270606;
        dwMax   = 1343.39;
        [~, rfmSignal, rfpSignal, rffSignal, rfBW] = ...
            sequence.waveforms.adiabaticHypSec(rfPulse.duration,...
            b1Max,dwMax,dt,gamma);
        
    case 'tantanh'
        rfPulse.doSliceSelect = 0;
        % predefined values
        b1Max   = 0.000015;
        dwMax   = 9500;
        beta    = 10;
        kappa   = 1.525373047373320;
        [~, rfmSignal, rfpSignal, rffSignal, rfBW] = ...
            sequence.waveforms.adiabaticTanTanhFull(rfPulse.duration,...
            b1Max,dwMax,beta,kappa,dt,gamma);
        
    case 'bir4'
        rfPulse.doSliceSelect = 0;
        % predefined values
        b1Max   = 0.000015;
        dwMax   = 8489;
        beta    = 5;
        kappa   = atan(10); %1.4711;
        [~, rfmSignal, rfpSignal, rffSignal, rfBW] = ...
            sequence.waveforms.adiabaticBIR4(rfPulse.duration,...
            rfPulse.flipAngle,b1Max,dwMax,beta,kappa,dt,gamma);
    
    case 'gauss'
        [~, rfmSignal, rfpSignal, rffSignal, rfBW] = ...
            sequence.waveforms.rfGauss( rfPulse.duration,...
            rfPulse.flipAngle,[],dt,gamma);
        
    otherwise % sinc as default
        [~, rfmSignal, rfpSignal, rffSignal, rfBW] = ...
            sequence.waveforms.rfSinc( rfPulse.duration,...
            rfPulse.flipAngle,rfPulse.cycles,dt,gamma);

end
        
% slice selection
if rfPulse.doSliceSelect

    % STEP 2 - Design slice selection gradient
    %   compute the required amplitude
    grMagnitude = rfBW/(gamma*rfPulse.sliceThickness);
    
    % check if maximum gradient limit is trespassed in Slice selection
    if grMagnitude > maxGStrenght
        msg = sprintf( ['Gradient amplitude in readout (%.3fmT/m) '...
            'exceeds system maximum (%.3fmT/m). '...
            'Consider increasing the slice thickness.'],...
            grMagnitude*1e3, maxGStrenght*1e3 );
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    end

    %   and make sure plateau time is enough for rfTime
    plateauTime = length(rfmSignal)*dt;
    [~,grSignal,grArea,~,rfLimits] = sequence.waveforms.grTrapPlateau(...
        plateauTime,grMagnitude,gradSlewrate,dt);

    % total number of time points
    numGrPoints = length(grSignal);
    numTimePoints = numGrPoints;
    
    % PRE-Rewind gradient
    if ( abs(rfPulse.preRWScale) > 1e-3 ) ...
            || ( abs(rfPulse.postRWScale) > 1e-3 ) ...
            || ( rfPulse.keepTimingSS > 0 )
        
        % generate RW gradient with half the area
        [~,rwSignal,~,~,~] = sequence.waveforms.grTrapArea(grArea/2,...
            0.99*maxGStrenght,gradSlewrate,dt);
        
        numRwPoints = length(rwSignal);
        
        % correct rf limits and the number of time points
        if ( abs(rfPulse.preRWScale) > 1e-3 ) || (rfPulse.keepTimingSS > 0)
            numTimePoints = numTimePoints + numRwPoints;
            rfLimits = rfLimits + numRwPoints;
        end
        if ( abs(rfPulse.postRWScale) > 1e-3 ) || (rfPulse.keepTimingSS > 0)
            numTimePoints = numTimePoints + numRwPoints;
        end
        
    end
    
    %  allocate signals
    rfm     = zeros(numTimePoints,1);
    rfp     = zeros(numTimePoints,1);
    rff     = zeros(numTimePoints,1);
    gr      = zeros(numTimePoints,1);
    time    = zeros(numTimePoints,1);
    
    % assign
    time(:) = dt*(1:numTimePoints);
    % incorporate RF in correct position
    assert(diff(rfLimits)+1==length(rfmSignal),...
        sprintf('\nERROR: %s : gradient plateu and RF dimension missmatch \n',functionName))
    rfm(rfLimits(1):rfLimits(2)) = rfmSignal(:);
    rfp(rfLimits(1):rfLimits(2)) = rfpSignal(:) + rfPulse.phase;
    rff(rfLimits(1):rfLimits(2)) = rffSignal(:);
    % rewind and gradinet
    if ( abs(rfPulse.preRWScale) > 1e-3 ) || (rfPulse.keepTimingSS > 0)
        gr(1:numRwPoints) = -rfPulse.preRWScale*rwSignal(:);
        gr(numRwPoints+1:numRwPoints+numGrPoints) = grSignal(:);
        counter = numRwPoints+numGrPoints;
    else
        gr(1:numGrPoints) = grSignal(:);
        counter = numGrPoints;
    end
    % rewind
    if ( abs(rfPulse.postRWScale) > 1e-3 ) || (rfPulse.keepTimingSS > 0)
        gr(counter+1:numTimePoints) = -rfPulse.postRWScale*rwSignal(:);
    end
    
else
    
    % No slice selection, signals is just RF
    grMagnitude = 0;
    numTimePoints = length(rfmSignal);
    % start and end of RF
    rfLimits = [1,numTimePoints];
    % allocate   
    rfm     = zeros(numTimePoints,1);
    rfp     = zeros(numTimePoints,1);
    rff     = zeros(numTimePoints,1);
    gr      = zeros(numTimePoints,1);
    time    = zeros(numTimePoints,1);
    % assign
    time(:)                      = dt*(1:numTimePoints);
    rfm(rfLimits(1):rfLimits(2)) = rfmSignal(:);
    rfp(rfLimits(1):rfLimits(2)) = rfpSignal(:) + rfPulse.phase;
    rff(rfLimits(1):rfLimits(2)) = rffSignal(:);
    
end

% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fms',...
        functionName, toc(tTotal) );
    fprintf(fid, '\n  IRL Time         %.3fms', time(end)*1e3);
    fprintf(fid, '\n  RF  Time         %.3fms', rfPulse.duration*1e3);
    fprintf(fid, '\n  RF  Flip Angle   %.1f', rfPulse.flipAngle);
    fprintf(fid, '\n  RF  Phase        %.1f', rfPulse.phase);
    fprintf(fid, '\n  RF  BandWidth    %.3fMHz', rfBW*1e-3);
    fprintf(fid, '\n  Slice Thickness  %.1fmm', rfPulse.sliceThickness*1e3);
	fprintf(fid, '\n  Gradient Level   %.1fmT/m', grMagnitude*1e3);
    fprintf(fid, '\n  Peak Gradient    %.1fmT/m', max(abs(gr))*1e3);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
