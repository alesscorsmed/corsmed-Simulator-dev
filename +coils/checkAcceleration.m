function [needCoilConfirm] = checkAcceleration(coilSystem, acquisition, expControl)
% COILS.CHECKACCELERATION
%
%   Verifies that for SENSE acceleration, 
%   the number of coil elements in the selected coil
%   allows for a correct reconstruction.
%   If not, warns the user, but let's the experiment run.
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils:checkAcceleration';
if (nargin < 2)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 3) || isempty(expControl)
    expControl.debug.debugMode  = 0;
    expControl.debug.debugFile  = [];
    expControl.connLocalDB      = [];
end
%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

% output
needCoilConfirm  = 0;

% check if SENSE is applied
if ~strcmpi(acquisition.data.parallelImaging, 'no')
    rFactor = acquisition.data.rFactor;
else
    rFactor = 0;
end

% in the case of SENSE acceleration
if rFactor
    
    % get the name and index of the RX coil
    activeRxName    = coilSystem.activeRx;
    activeRxIndex   = coilSystem.indexRx;
    activeRxData    = coilSystem.coilModel{activeRxIndex}.data;
    activeRxNumCoil = activeRxData.numCoils;
    % extract acceleration capabilities of selected coil
    parallelAP = activeRxData.parallelAP; % max parallelism in AP direction (0 is no parallel)
    parallelRL = activeRxData.parallelRL; % max parallelism in RL direction (0 is no parallel)
    parallelFH = activeRxData.parallelFH; % max parallelism in FH direction (0 is no parallel)
    
    % get main plane
    %% find plane orientation
    pointsAll           = acquisition.data.pointsAll;
    pointsAllFrontend   = acquisition.data.pointsAllFrontend;
    % extract coordinates data
    points          = pointsAll{1};
    points          = strsplit(points,';');
    pointsFrontEnd  = pointsAllFrontend{1};
    pointsFrontEnd  = strsplit(pointsFrontEnd,';');
    % generate plane
    % assign points
    plane.LTop = str2num(points{1,1});  % Top-Left
    plane.RTop = str2num(points{1,2});  % Top-Right
    plane.LBot = str2num(points{1,3});  % Bottom-Left
    plane.RBot = str2num(points{1,4});  % Bottom-right
    % get the points from the Front End
    plane.LTopFrontEnd  = str2num(pointsFrontEnd{1,1});  % Top-Left
    plane.RTopFrontEnd  = str2num(pointsFrontEnd{1,2});  % Top-Right
    plane.LBotFrontEnd  = str2num(pointsFrontEnd{1,3});  % Bottom-Left
    plane.RBotFrontEnd  = str2num(pointsFrontEnd{1,4});  % Bottom-right
    % compute the transformations and final plane coordinates
    [plane] = domain.planeHandling.generatePlane(plane);

    % verify acceleration does not cause issues
    accPlane = sort([plane.BOrient, plane.TOrient]);
    accWarning = [];
    switch lower(accPlane)
        case 'ap'
            if rFactor >= parallelAP
                accWarning = sprintf(['Acceleration factor R=%d ',...
                    'may cause artefacts for Coil %s with %d elements.'],...
                    rFactor, activeRxName, activeRxNumCoil );
            end
        case 'lr'
            if rFactor >= parallelRL
                accWarning = sprintf(['Acceleration factor R=%d ',...
                    'may cause artefacts for Coil %s with %d elements.'],...
                    rFactor, activeRxName, activeRxNumCoil );
            end
        case 'fh'
            if rFactor >= parallelFH
                if (rFactor < parallelRL) || (rFactor < parallelAP)
                    accWarning = sprintf(['Acceleration factor R=%d in FH direction',...
                        'may cause artefacts for Coil %s with %d elements. '...
                        'You may considering changing the encoding direction.'],...
                        rFactor, activeRxName, activeRxNumCoil );
                else
                    accWarning = sprintf(['Acceleration factor R=%d ',...
                        'may cause artefacts for Coil %s with %d elements.' ],...
                        rFactor, activeRxName, activeRxNumCoil );
                end
            end
        otherwise
            accWarning = [];
    end
    
    if ~isempty(accWarning)
        needCoilConfirm = 1;
        if isfield(expControl,'application') && strcmpi(expControl.application, 'edutool')
            % EduTool warning
            eduTool.frontend.updateExperimentProgress(expControl,'','confirm',accWarning);
        end
        % Screen warning
        fprintf(1, '\n');
        fprintf(1,'WARNING: %s\n', accWarning);
        fprintf(1, '\n');
    end
    
end


%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done',functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Active coil       %s', coilSystem.activeRx);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end