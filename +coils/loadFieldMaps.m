function [coilSystem] = loadFieldMaps(coilSystem)
%
% COILS.LOADFIELDMAPS
%
%     loads the active coil field maps
%
% INPUT
%   coilSystem          struct with coils data
%
% OUTPUT
%   coilSystem          updated struct with correct active coils
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.loadFieldMaps';
if (nargin < 1)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

try 
    %% load the required field data for Tx coil
    indexTx = coilSystem.indexTx;
    if isempty(coilSystem.coilModel{indexTx}.maps)
        coilSystem.coilModel{indexTx}.maps = load( ...
            coilSystem.coilModel{indexTx}.data.mapsFile,...
            'bx', 'by', 'sosB1m', 'sosB1p', 'sosErms', 'dim', 'spatial' );
    end    
catch
    %% no data available, error
    message = sprintf('WARNING at %s: no EM fields data available for %s TX coil',...
        functionName, coilSystem.coilModel{indexTx}.data.name );
    fprintf(1, '\n\n%s\n\n', message);
    %% switch to default
    [coilSystem] = coils.updateActiveCoils(coilSystem.defaultTx, coilSystem);
end

try 
    %% load the required field data for Rx coil
    indexRx = coilSystem.indexRx;
    if isempty(coilSystem.coilModel{indexRx}.maps)
        coilSystem.coilModel{indexRx}.maps = load( ...
            coilSystem.coilModel{indexRx}.data.mapsFile,...
            'bx', 'by', 'sosB1m', 'sosB1p', 'sosErms', 'dim', 'spatial' );
    end 
catch
    %% no data available, error
    message = sprintf('WARNING at %s: no EM fields data available for %s RX coil',...
        functionName, coilSystem.coilModel{indexRx}.data.name );
    fprintf(1, '\n\n%s\n\n', message);
    %% switch to default
    [coilSystem] = coils.updateActiveCoils(coilSystem.defaultTx, coilSystem);
end


