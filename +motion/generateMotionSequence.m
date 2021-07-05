function [motion] = generateMotionSequence(sequence,numNEX,expControl)
%
% MOTION.GENERATEMOTIONSEQUENCE
%
%	Generates a motion sequence from motion specs.
%
% INPUT
%   sequence        sequence with time signals
%   expControl      experiment control struct
%
% OUTPUT
%   motion        struct with time domain motion transformations
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'motion.generateMotionSequence';
if (nargin < 2)
    ME = MException('simulator:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% initialize motion
motion.type = 'none';
motion.angleRotXY   = 0.0;
motion.xCenterRotXY = 0.0;
motion.yCenterRotXY = 0.0;
% assign ZERO rotation in XZ
motion.angleRotXZ   = 0.0;
motion.xCenterRotXZ = 0.0;
motion.zCenterRotXZ = 0.0;
% assign ZERO rotation in YZ
motion.angleRotYZ   = 0.0;
motion.yCenterRotYZ = 0.0;
motion.zCenterRotYZ = 0.0;
% assign ZERO translation
motion.transX = 0.0;
motion.transY = 0.0;
motion.transZ = 0.0;

%% create the motion depending on the pattern
motionSpecs = expControl.motionSpecs;

% current time shift for the Sequence NEX number
timeShift = (numNEX-1)*sequence.totalTime;

switch lower(motionSpecs.pattern)
    case 'rotational'
        %% periodic rotation on XY plane
        % create a time dependent angle between -rotAngle and +rotAngle ??
        % varying with frequency rotFreq ??
        time        = sequence.time(:) + timeShift; % sequence time
        rotAngle    = motionSpecs.rotAngle*pi/180; % in radians
        rotFreq     = motionSpecs.rotFreq;
        % periodic movement between 0 and 2
        % angleArray = 1 + sin(2*pi*rotFreq*time); % sinusoidal
        angleArray = 1 + sawtooth(2*pi*rotFreq*time,1/2); % linear
        
        % assign rotation in XY
        motion.type = 'XY rotation';
        % periodic between 0 and rotAngle
        motion.angleRotXY   = 0.5*rotAngle*angleArray;
        motion.xCenterRotXY = 0.0;
        motion.yCenterRotXY = 0.0;
        % assign ZERO rotation in XZ
        motion.angleRotXZ   = 0*angleArray;
        motion.xCenterRotXZ = 0.0;
        motion.zCenterRotXZ = 0.0;
        % assign ZERO rotation in YZ
        motion.angleRotYZ   = 0*angleArray;
        motion.yCenterRotYZ = 0.0;
        motion.zCenterRotYZ = 0.0;
        
        % assign ZERO translation
        motion.transX = 0*angleArray;
        motion.transY = 0*angleArray;
        motion.transZ = 0*angleArray;
        
    case 'translational'
        %% periodic translation in one of the axis
        % create a time dependent translation
        % varying with frequency rotFreq ??
        time        = sequence.time(:) + timeShift; % sequence time
        transMag    = motionSpecs.transMag;
        transFreq   = motionSpecs.transFreq;
        % periodic movement between 0 and 2
        % transArray = 1 + sin(2*pi*transFreq*time); % sinusoidal
        transArray = 1 + sawtooth(2*pi*transFreq*time,1/2); % linear
        
        % assign ZERO rotation in XY
        motion.angleRotXY   = 0*transArray;
        motion.xCenterRotXY = 0.0;
        motion.yCenterRotXY = 0.0;
        % assign ZERO rotation in XZ
        motion.angleRotXZ   = 0*transArray;
        motion.xCenterRotXZ = 0.0;
        motion.zCenterRotXZ = 0.0;
        % assign ZERO rotation in YZ
        motion.angleRotYZ   = 0*transArray;
        motion.yCenterRotYZ = 0.0;
        motion.zCenterRotYZ = 0.0;
        
        % translation
        switch motionSpecs.transAxis            
            case 1 % X
                motion.type = 'X translation';
                % assign ZERO translation
                motion.transX = 0.5*transMag*transArray;
                motion.transY = 0*transArray;
                motion.transZ = 0*transArray;
            case 2 % Y
                motion.type = 'Y translation';
                % assign ZERO translation
                motion.transX = 0*transArray;
                motion.transY = 0.5*transMag*transArray;
                motion.transZ = 0*transArray;    
            case 3 % Z
                motion.type = 'Z translation';
                % assign ZERO translation
                motion.transX = 0*transArray;
                motion.transY = 0*transArray;
                motion.transZ = 0.5*transMag*transArray;   
            otherwise % no movement
                motion.type = 'none';      
        end
        
    otherwise
        %% No motion
        motion.type = 'none';   
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done, elapsed time %.3fs',...
        functionName, tTotal);
    fprintf(fid, '\n  Motion Type:  %s', motion.type);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
