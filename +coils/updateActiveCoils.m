function [coilSystem] = updateActiveCoils(coilName, coilSystem)
%
% COILS.UPDATEACTIVECOILS
%
%     Updates the active Tx and Rx coils to match the given coilName.
%
% INPUT
%   coilName        name of the coil to activate
%   coilSystem      struct with coils data
%
% OUTPUT
%   coilSystem      updated struct with correct active coils
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.updateActiveCoils';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% if desired coil is not currently active, update active coils
if ~strcmpi(coilName, coilSystem.activeRx)
    coilSystem.indexRx = 0; % set to zero to try to find coil
    for cc = 1:coilSystem.numModels
        if strcmpi(coilName, coilSystem.coilModel{cc}.data.name)
            coilSystem.activeRx     = coilSystem.coilModel{cc}.data.name;
            coilSystem.indexRx      = cc;
            if (coilSystem.coilModel{cc}.data.isTx)
                coilSystem.activeTx     = coilSystem.coilModel{cc}.data.name;
                coilSystem.indexTx      = cc;
            else
                coilSystem.activeTx     = coilSystem.defaultTx;
                coilSystem.indexTx      = coilSystem.defaultIndex;
            end
            break;
        end
    end
end

%% check that a valid coil was found
if (coilSystem.indexRx == 0)
    % use the default otherwise, and report
    coilSystem.activeRx     = coilSystem.defaultTx;
    coilSystem.indexRx      = coilSystem.defaultIndex;
    coilSystem.activeTx     = coilSystem.defaultTx;
    coilSystem.indexTx      = coilSystem.defaultIndex;
    message = sprintf('WARNING at %s: coil %s not found, default (%s) will be used',...
        functionName, coilName, coilSystem.defaultTx );
    fprintf(1, '\n\n%s\n\n', message);
end