function [pulseSequence] = generatePulseSequence( ...
         acquisition, encoding, mrSystem, expControl, anatomicalModel )
%
% SEQUENCE.GENERATEPULSESEQUENCE
%
%	Interface to generate pulse sequences based on.
%
% INPUT
%   acquisition         
%   encoding
%   mrSystem        
%   expControl      
%
% OUTPUT
%   pulseSequence   pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.generatePulseSequence';

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

%% generate the required sequence
switch lower(acquisition.data.pulseSeqFamilyName)
    
    case {lower('GRE'),lower('perfusion-gre')}
        
        % modify RF pulse to comply with case: unbalanced refocusing
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % define no spoiling:
        %% CHANGE: Correct
        acquisition.data.spoiled =1;
        
        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE2D(...
            acquisition, encoding, mrSystem, expControl);
        
    case lower('SWI-2D-GRE')
        
        % modify RF pulse to comply with case: unbalanced refocusing
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % define no spoiling:
        %% CHANGE: Correct
        acquisition.data.spoiled =1;
        
        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE2D(...
            acquisition, encoding, mrSystem, expControl);
        
    case lower('spoiled-GRE') 
        
        % modify RF pulse to comply with case: unbalanced refocusing
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % define spoiling:
        acquisition.data.spoiled = 1;

        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE2D(...
            acquisition, encoding, mrSystem, expControl);
        
    case lower('GRE-3D')
        
        % no slice selection
        acquisition.mainRF.doSliceSelect  = 0;
        
        % define no spoiling:
        acquisition.data.spoiled = 1;
        
        % deactivate slice selection for 3D
        expControl.sequence.deactivateSS = 1;
        
        % update Slice (3D) ecoding data
        acquisition.data.fovSE = acquisition.data.numSlices ...
            * acquisition.data.sliceThickness;
        acquisition.data.numSE = acquisition.data.numSlices;
        
        % will use the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE3D(...
            acquisition, encoding, mrSystem, expControl);

    case lower('spoiled-GRE-3D')
        
        % no slice selection
        acquisition.mainRF.doSliceSelect  = 0;
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % deactivate slice selection for 3D
        expControl.sequence.deactivateSS = 1;
 
        % update Slice (3D) ecoding data
        acquisition.data.fovSE = acquisition.data.numSlices ...
            * acquisition.data.sliceThickness;
        acquisition.data.numSE = acquisition.data.numSlices;
        
        % will use the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE3D(...
            acquisition, encoding, mrSystem, expControl);
        
    case lower('GRE-OOP')
        
        % modify RF pulse to comply with case: unbalanced refocusing
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % define no spoiling:
        %% CHANGE: Correct
        acquisition.data.spoiled =1;
        
        % will ignore the Slice (3D) encodings
        [pulseSequence] = sequence.familyGRE.GRE2D(...
            acquisition, encoding, mrSystem, expControl);
    
    case lower('MP-RAGE')
        
        % update Slice (3D) ecoding data
        acquisition.data.fovSE = acquisition.data.numSlices ...
            * acquisition.data.sliceThickness;
        acquisition.data.numSE = acquisition.data.numSlices;
        
        % no slice selection
        acquisition.mainRF.doSliceSelect  = 0;
        
        % force TR of each encoding to shortest possible
        acquisition.data.forceMinTR = 0; % no, TR given by user
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % deactivate slice selection for 3D
        expControl.sequence.deactivateSS = 1;
        
        % IR-prep
        acquisition.prepIR.Apply = 1; % apply IR
        
        % call function
        [pulseSequence] = sequence.familyGRE.MPRAGE(...
            acquisition, encoding, mrSystem, expControl);
    
    case lower('GRE-conc')
        
        % same as 2D
        % modify RF pulse to comply with case: unbalanced refocusing
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % force TR of each encoding to shortest possible
        acquisition.data.forceMinTR = 1; % to concatenate
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % do not deactivate slice selection: multi-slice for 3D
        expControl.sequence.deactivateSS = 0;
        
        % call function
        [pulseSequence] = sequence.familyGRE.sliceHopGRE(...
            acquisition, encoding, mrSystem, expControl);        

    case {lower('bSSFP'),lower('cine-bssfp')}
        
        % define RF properties
        % keepTimingSS: 
        %   keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        if strcmpi(acquisition.data.fatsat, 'lipidir')
            % call the IR function
            [pulseSequence] = sequence.familySSFP.irBalancedSSFP(...
                acquisition, encoding, mrSystem, expControl);
        else
            % call the function
            [pulseSequence] = sequence.familySSFP.balancedSSFP(...
                acquisition, encoding, mrSystem, expControl);
        end
        
    case lower('IR-bSSFP')
        
        % define RF properties
        % keepTimingSS: 
        %   keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % call the function
        [pulseSequence] = sequence.familySSFP.irBalancedSSFP(...
            acquisition, encoding, mrSystem, expControl);
    
    case lower('MOLLI')
        
        % define RF properties
        % keepTimingSS: 
        %   keeps the timing of the pre and post RW, even if 0,
        %   to maintaing symmetry
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % call the function
        [pulseSequence] = sequence.familySSFP.MOLLI(...
            acquisition, encoding, mrSystem, expControl, anatomicalModel);
    
        
    case lower('SE') % 'SpinEcho2D'
        
        % main RF
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        % acquisition.refRF.flipAngle      = 180; % comes from DB
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % call function
        [pulseSequence] = sequence.familySE.spinEcho2D(...
            acquisition, encoding, mrSystem, expControl);
    
    case lower('IR-SE') % 'IR-SpinEcho2D'
        
        % main RF
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        % acquisition.refRF.flipAngle      = 180; % comes from DB
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % apply IR
        acquisition.prepIR.Apply = 1; % apply IR
        
        % call general spin echo function
        [pulseSequence] = sequence.familySE.spinEcho(...
            acquisition, encoding, mrSystem, expControl);
        
     case lower('SE-3D') % 'SpinEcho3D'
        
        % main RF, no slice selection (3D)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
 
        % update Slice (3D) ecoding data
        acquisition.data.fovSE = acquisition.data.numSlices ...
            * acquisition.data.sliceThickness;
        acquisition.data.numSE = acquisition.data.numSlices;      
        
        % call function
        [pulseSequence] = sequence.familySE.spinEcho3D(...
            acquisition, encoding, mrSystem, expControl);    
        
    case lower('IR-SE-3D') % 'SpinEcho3D'
        
        % main RF, no slice selection (3D)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % apply IR
        acquisition.prepIR.Apply = 1; % apply IR
 
        % update Slice (3D) ecoding data
        acquisition.data.fovSE = acquisition.data.numSlices ...
            * acquisition.data.sliceThickness;
        acquisition.data.numSE = acquisition.data.numSlices;
                
        % call general spin echo function
        [pulseSequence] = sequence.familySE.spinEcho(...
            acquisition, encoding, mrSystem, expControl);

    case lower('SS-FSE')
        
        % main RF
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        % acquisition.refRF.flipAngle      = 180; % comes from DB
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        if strcmpi(acquisition.data.fatsat, 'lipidir')
            % apply the IR
            acquisition.prepIR.Apply = 1;
        else
            % do not apply IR
            acquisition.prepIR.Apply = 0;
        end
        
        % call function
        [pulseSequence] = sequence.familySE.singleShotFSE(...
            acquisition, encoding, mrSystem, expControl);

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
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % apply IR
        acquisition.prepIR.Apply = 1; % apply IR
        
        % call function
        [pulseSequence] = sequence.familySE.singleShotFSE(...
            acquisition, encoding, mrSystem, expControl);

    case lower('TSE')
        
        % main RF
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % Refocusing RF: 180
        % acquisition.refRF.flipAngle      = 180; % comes from DB
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 1;
        acquisition.refRF.preRWScale     = 1; % rewind SS before RF
        acquisition.refRF.postRWScale    = 1; % rewind SS after RF
        acquisition.refRF.keepTimingSS   = 1;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        if strcmpi(acquisition.data.fatsat, 'lipidir')
            % apply the IR
            acquisition.prepIR.Apply = 1;
        else
            % do not apply IR
            acquisition.prepIR.Apply = 0;
        end
        
        % call function
        [pulseSequence] = sequence.familySE.turboSE(...
            acquisition, encoding, mrSystem, expControl);

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
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % apply IR
        acquisition.prepIR.Apply = 1; % apply IR
        
        % call function
        [pulseSequence] = sequence.familySE.turboSE(...
            acquisition, encoding, mrSystem, expControl);

    case lower('EPI')  
        
        acquisition.mainRF.doSliceSelect  = 1; % slice selection
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % call function
        [pulseSequence] = sequence.familyEPI.greEPI(...
            acquisition, encoding, mrSystem, expControl);
        
    case lower('SE-EPI')  
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % call function
        [pulseSequence] = sequence.familyEPI.spinEchoEPI(...
            acquisition, encoding, mrSystem, expControl);
        
    case lower('PG-SE-EPI')  
        
        % main RF
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = -pi/2;
        acquisition.mainRF.doSliceSelect  = 0;

        % Refocusing RF: 180
        acquisition.refRF.flipAngle      = 180;
        acquisition.refRF.phase          = 0;
        acquisition.refRF.doSliceSelect  = 0;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 1;
        
        % call function
        [pulseSequence] = sequence.familyPGSE.pulsedGradientSpinEchoEPI(...
            acquisition, encoding, mrSystem, expControl);
        
        % assign updated info about PGSE
        acquisition.encPG = pulseSequence.encPG;
    
    
    case lower('bssfp-train')  
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 90;
        acquisition.mainRF.phase          = 0;
        
        acquisition.mainRF.doSliceSelect  = 1;
        acquisition.mainRF.preRWScale     = 1; % rewind SS before RF
        acquisition.mainRF.postRWScale    = 1; % rewind SS after RF
        acquisition.mainRF.keepTimingSS   = 1;
        
        % deactivate slice selection
        expControl.sequence.deactivateSS = 0;
        
        % IR-prep
        acquisition.prepIR.Apply = 1; % apply IR
        
        % use a cos-modulated train of pulses (arbitrary)
        expControl.sequence.rfCosTrain   = 1;
        
        % skip preparation after inversion (not SS anyway)
        acquisition.data.dummySSFP       = 0;
        
        % call the function
        [pulseSequence] = sequence.familySSFP.irBalancedSSFP(...
            acquisition, encoding, mrSystem, expControl);
    
    
    case lower('fisp-train')  
        
        % RF type (more parametes available)
        acquisition.mainRF.flipAngle      = 45;
        acquisition.mainRF.phase          = 0;
        
        acquisition.mainRF.preRWScale     = 0; % no rewind before RF (spoiler)
        acquisition.mainRF.postRWScale    = 1; % rewind after RF
        acquisition.mainRF.keepTimingSS   = 0; % do not keep timings
        
        % force TR of each encoding to shortest possible
        acquisition.data.forceMinTR = 1; % to concatenate
        
        % define spoiling:
        acquisition.data.spoiled = 1;
        
        % do not deactivate slice selection: multi-slice for 3D
        expControl.sequence.deactivateSS = 0;
        
        % deactivate slice selection
        acquisition.mainRF.doSliceSelect  = 1;
        
        % skip preparation after inversion (not SS anyway)
        acquisition.data.dummySSFP       = 0;
        
        % use a cos-modulated train of pulses (arbitrary)
        expControl.sequence.rfCosTrain   = 1;
        
        % build it
        [pulseSequence] = sequence.familyGRE.irFISP(...
            acquisition, encoding, mrSystem, expControl);
        
        
    otherwise
        
        msg = sprintf( ['The selected sequence type (%s) is not available. ',...
            'If error persist please contact the Corsmed team.'],...
            acquisition.data.pulseSeqFamilyName );
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
        
end


%% apply minimization of kernel context change 
%  only if there was time compression
%  this reduces the number of parts of the sequence 
%  by combining compatible parts.
if expControl.sequence.timeCompression
    [pulseSequence] = sequence.tools.minContextChange(...
        pulseSequence, expControl);
end

%% compute time differences
pulseSequence.timeDiff = [pulseSequence.time(1); diff(pulseSequence.time(:))];

%% family name and number, and NEX
pulseSequence.seqNum     			= acquisition.data.pulseSeqNum;
pulseSequence.familyName 			= acquisition.data.pulseSeqFamilyName;
pulseSequence.NEX        			= acquisition.data.NEX;
pulseSequence.pulseSeqActualName    = acquisition.data.pulseSeqActualName;
if strcmp(pulseSequence.pulseSeqActualName,'DIXON')
    pulseSequence.displaySeqName   	= pulseSequence.pulseSeqActualName;
else
    pulseSequence.displaySeqName   	= pulseSequence.familyName;
end

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
