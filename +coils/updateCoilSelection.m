function [coilSystem] = updateCoilSelection(expControl,anatomicalModel,coilSystem)
% COILS.UPDATECOILSELECTION
%
%     Updates the selected coil, positioning and SAR values.
%
% INPUT
%   expControl      control struct with experiment info
%   anatomicalModel struct with model data
%   coilSystem      struct with coils data
%
% OUTPUT
%   coilSystem      updated struct with correct active coils
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.updateCoilSelection';
if (nargin < 3)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
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
%% extract relevant info
coilName   = expControl.model.coilType;
zIsocenter = expControl.model.zIsocenter;
%% update the active coils, both Tx and Rx
[coilSystem] = coils.updateActiveCoils(coilName, coilSystem);
%% load field maps if not there
[coilSystem] = coils.loadFieldMaps(coilSystem);
%% check if need to perform update in SAR
doUpdate = isempty( coilSystem.coilModel{coilSystem.indexTx}.sar )...
    || ( abs(coilSystem.isocenter(3) - zIsocenter) > 1e-4 );
if doUpdate
    %% update isocenter
    coilSystem.isocenter(3) = zIsocenter;
    %% compute SAR for the transmit coil with the defined isocenter
    [coilSystem] = coils.computeSAR(coilSystem, anatomicalModel);
end
%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Active coil       %s', coilSystem.activeRx);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end