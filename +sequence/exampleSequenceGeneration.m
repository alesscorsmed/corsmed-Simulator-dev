% 
% Example of how to generate sequences
%
%========================  CORSMED AB © 2020 ==============================
%

clear all
close all
clc

acquisition     = data.acquisition.initialize();
expControl      = data.expControl.initialize();
mrSystem        = data.mrSystem.initialize();

mrSystem.maxGStrenght   = 0.040; % controls maximum gradient
mrSystem.SlewRate       = 150.0; % set to 0 for ideal grads with no ramps

acquisition.data.fovFE       = 0.2; % field of view in m
acquisition.data.fovPE       = 0.2; % field of view in m
acquisition.data.fovSE       = 0.2; % field of view in m
acquisition.data.numFE       = 256; % freq encodings
acquisition.data.numPE       = 10; % phase encodings
acquisition.data.numSE       = 5; % 3D encodings
acquisition.data.TR          = 2000e-3;
acquisition.data.TE          = 5e-3;
acquisition.prepIR.TI        = 1000e-3;
% acceleration
acquisition.data.parallelImaging    = 'no';
acquisition.data.rFactor            = 2;
acquisition.data.partialFourier     = 'none';%'phaseConjugate';
acquisition.data.fFactor            = 1.0;
% used only in SS-FSE, force TE to shortest possible
acquisition.data.forceMinTE  = 0;
acquisition.mainRF.duration = 1e-3;
acquisition.refRF.duration  = 1e-3;
% used in MP-RAGE: force TR to shortest possible
acquisition.data.forceMinTR = 0;
% ETL for multi-shot
acquisition.data.ETL = 1;

acquisition.data.pulseSeqFamilyName = 'MP-RAGE';
acquisition.data.is3D = 1;
acquisition.data.numSlices = acquisition.data.numSE;

[kSpaceInfo] = sequence.tools.cartesianEncodingInfo(...
    acquisition.data, expControl);
acquisition.kSpaceInfo = kSpaceInfo;

tTotal = tic();

switch lower(acquisition.data.pulseSeqFamilyName)
    
    case lower('GRE2D') 
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 15;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 0;
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        % keepTimingSS
        %   if 1 keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.keepTimingSS   = 0;
        
        % define no spoiling:
        acquisition.data.spoiled = 0;
        
        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE2D(...
            acquisition, mrSystem, expControl);
    
    case lower('SpoiledGRE2D') 
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 15;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 0;
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        % keepTimingSS
        %   if 1 keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.keepTimingSS   = 0;
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE2D(...
            acquisition, mrSystem, expControl);
        
    case lower('GRE3D')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 30;
        acquisition.mainRF.doSliceSelect  = 0;
        
        % define no spoiling:
        acquisition.data.spoiled = 0;
        
        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE3D(...
            acquisition, mrSystem, expControl);
        
    case lower('SpoiledGRE3D')
        
        % no slice selection
        acquisition.mainRF.flipAngle      = 30;
        acquisition.mainRF.doSliceSelect  = 0;
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % will use the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE3D(...
            acquisition, mrSystem, expControl);
        
        
    case lower('MP-RAGE')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 15;
        acquisition.mainRF.doSliceSelect  = 0;
        acquisition.mainRF.preRWScale     = 0;
        acquisition.mainRF.postRWScale    = 0; % rewind after RF
        % keepTimingSS
        %   if 1 keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.keepTimingSS   = 0;
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % force TR of each encoding to shortest possible
        acquisition.data.forceMinTR = 1;
        
        % IR-prep
        acquisition.prepIR.Apply = 1; % apply IR
        
        % call function
        [pulseSequence] = sequence.familyGRE.MPRAGE(...
            acquisition, mrSystem, expControl);
        
        
    case lower('bSSFP')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 30;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        % keepTimingSS
        %   if 1 keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.keepTimingSS   = 1;
        
        acquisition.data.dummySSFP = 4; % dummy reps as preparation for steady-state
        
        % will ignore the SE encodings
        [pulseSequence] = sequence.familySSFP.balancedSSFP(...
            acquisition, mrSystem, expControl);
    
    case lower('IR-bSSFP')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 30;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        % keepTimingSS
        %   if 1 keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.keepTimingSS   = 1;
        
        acquisition.data.dummySSFP = 4; % dummy reps as preparation for steady-state
        
        % will ignore the SE encodings
        [pulseSequence] = sequence.familySSFP.irBalancedSSFP(...
            acquisition, mrSystem, expControl);
        
    case lower('SpinEcho2D')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % call function
        [pulseSequence] = sequence.familySE.spinEcho2D(...
            acquisition, mrSystem, expControl);
        
    case lower('SpinEcho3D')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % call function
        [pulseSequence] = sequence.familySE.spinEcho3D(...
            acquisition, mrSystem, expControl);

    case lower('SS-FSE')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % call function
        [pulseSequence] = sequence.familySE.singleShotFSE(...
            acquisition, mrSystem, expControl);

    case lower('IR-SS-FSE')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % IR
        acquisition.prepIR.Apply = 1; % do not apply IR
        
        % call function
        [pulseSequence] = sequence.familySE.singleShotFSE(...
            acquisition, mrSystem, expControl);

        
    case lower('TSE')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % No IR
        acquisition.prepIR.Apply = 0; % do not apply IR
        
        % call function
        [pulseSequence] = sequence.familySE.turboSE(...
            acquisition, mrSystem, expControl);
        
    case lower('IR-TSE')
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % IR
        acquisition.prepIR.Apply = 1; % do not apply IR
        
        % call function
        [pulseSequence] = sequence.familySE.turboSE(...
            acquisition, mrSystem, expControl);
        
    case lower('EPI')  
        
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % call function
        [pulseSequence] = sequence.familyEPI.greEPI(...
            acquisition, mrSystem, expControl);
    
    case lower('SE-EPI')  
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % call function
        [pulseSequence] = sequence.familyEPI.spinEchoEPI(...
            acquisition, mrSystem, expControl);
        
    case lower('PG-SE-EPI')  
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % call function
        [pulseSequence] = sequence.familyPGSE.pulsedGradientSpinEchoEPI(...
            acquisition, mrSystem, expControl);
        
        % assign updated info about PGSE
        acquisition.encPG = pulseSequence.encPG;
        
        % apply the DW gradients
        pulseSequence.gxSignal = pulseSequence.gxSignal +...
            pulseSequence.gdwSignal(:,1);
        pulseSequence.gySignal = pulseSequence.gySignal +...
            pulseSequence.gdwSignal(:,2);
        pulseSequence.gzSignal = pulseSequence.gzSignal +...
            pulseSequence.gdwSignal(:,3);

    otherwise

end

fprintf(1,'\n\n%s sequence done, elapsed time %.3fs\n', ...
    acquisition.data.pulseSeqFamilyName, toc(tTotal));


%% plot
figure();
plot(pulseSequence.time,pulseSequence.rfmSignal*1e3);
xlabel('time (s)');
hold on
plot(pulseSequence.time,pulseSequence.gxSignal);
plot(pulseSequence.time,pulseSequence.gySignal);
plot(pulseSequence.time,pulseSequence.gzSignal);
plot(pulseSequence.time(pulseSequence.rxSignal>0), ...
    pulseSequence.gxSignal(pulseSequence.rxSignal>0), 'o');
plot(pulseSequence.time(pulseSequence.rxLimits(:,1)), ...
    pulseSequence.gxSignal(pulseSequence.rxLimits(:,1)), '>');
plot(pulseSequence.time(pulseSequence.rxLimits(:,2)), ...
    pulseSequence.gxSignal(pulseSequence.rxLimits(:,2)), '<');
plot(pulseSequence.time(pulseSequence.partLimits(:,1)),zeros(pulseSequence.numParts,1), '^');
plot(pulseSequence.time(pulseSequence.partLimits(:,2)),zeros(pulseSequence.numParts,1), 'v');
if nnz(pulseSequence.swcSignal) > 0
    plot(pulseSequence.time(pulseSequence.swcSignal>0), ...
        zeros(nnz(pulseSequence.swcSignal),1), 's', 'LineWidth', 2, 'MarkerSize', 10);
    legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end', 'SWC');
else
    legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end');
end

